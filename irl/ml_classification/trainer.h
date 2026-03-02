#pragma once

#include <torch/torch.h>
#include <vector>
#include <iostream>
#include <cstdint>

namespace IRL {

template<typename NetType, typename Loader>
class Trainer {
public:
    Trainer(NetType& net,
            torch::optim::Optimizer& optimizer,
            Loader* train_loader,
            Loader* test_loader,
            Loader* val_loader,
            int epochs,
            // Logging / outputs (owned by caller)
            std::vector<double>* train_loss,
            std::vector<double>* val_loss,
            std::vector<double>* test_accuracy_vec,
            std::vector<double>* val_accuracy_vec,
            double* test_loss_out,
            double* test_accuracy_out,
            std::vector<std::vector<int64_t>>* confusion_matrix)
        : net_(net),
          optimizer_(optimizer),
          train_loader_(train_loader),
          test_loader_(test_loader),
          val_loader_(val_loader),
          epochs_(epochs),
          train_loss_(train_loss),
          val_loss_(val_loss),
          test_accuracy_vec_(test_accuracy_vec),
          val_accuracy_vec_(val_accuracy_vec),
          test_loss_out_(test_loss_out),
          test_accuracy_out_(test_accuracy_out),
          confusion_matrix_(confusion_matrix) {}

    void train() {
        // Clear logs if provided
        if (train_loss_) train_loss_->clear();
        if (val_loss_) val_loss_->clear();
        if (test_accuracy_vec_) test_accuracy_vec_->clear();
        if (val_accuracy_vec_) val_accuracy_vec_->clear();
        if (test_loss_out_) *test_loss_out_ = 0.0;
        if (test_accuracy_out_) *test_accuracy_out_ = 0.0;

        for (int epoch = 1; epoch <= epochs_; ++epoch) {
            std::cout << "Epoch " << epoch << "/" << epochs_ << std::endl;

            net_.train();

            double loss_sum = 0.0;   // sum of per-sample losses over the epoch
            int64_t total_samples = 0;

            for (auto& batch : *train_loader_) {
                optimizer_.zero_grad();

                torch::Tensor prediction = net_.forward(batch.data);

                // cross_entropy returns mean over the batch by default
                torch::Tensor loss = torch::nn::functional::cross_entropy(prediction, batch.target);

                loss.backward();
                optimizer_.step();

                const int64_t bs = batch.target.size(0);
                loss_sum += loss.template item<double>() * static_cast<double>(bs);
                total_samples += bs;
            }

            const double mean_train_loss =
                (total_samples > 0) ? (loss_sum / static_cast<double>(total_samples)) : 0.0;

            if (train_loss_) train_loss_->push_back(mean_train_loss);

            // Validation: record loss + accuracy
            double mean_vloss = 0.0;
            double vacc = 0.0;
            evaluate(*val_loader_, &mean_vloss, &vacc);

            if (val_loss_) val_loss_->push_back(mean_vloss);
            if (val_accuracy_vec_) val_accuracy_vec_->push_back(vacc);

            // Test: record accuracy per epoch if requested (no printing)
            double mean_tloss = 0.0;
            double tacc = 0.0;
            evaluate(*test_loader_, &mean_tloss, &tacc);
            if (test_accuracy_vec_) test_accuracy_vec_->push_back(tacc);
        }

        // Final test evaluation + print
        double final_test_loss = 0.0;
        double final_test_acc = 0.0;
        evaluate(*test_loader_, &final_test_loss, &final_test_acc);

        if (test_loss_out_) *test_loss_out_ = final_test_loss;
        if (test_accuracy_out_) *test_accuracy_out_ = final_test_acc;

        std::cout << "Final Test Loss: " << final_test_loss
                  << " , Final Test Accuracy: " << final_test_acc << std::endl;

        if (confusion_matrix_) {
            // Infer number of classes from network output on one batch (safe fallback)
            int64_t num_classes = 0;
            {
                net_.eval();
                torch::NoGradGuard no_grad;
                for (auto& batch : *test_loader_) {
                    auto pred = net_.forward(batch.data);
                    num_classes = pred.size(1);
                    break;
                }
            }

            confusion_matrix_->assign(
                static_cast<size_t>(num_classes),
                std::vector<int64_t>(static_cast<size_t>(num_classes), 0)
            );

            // Fill confusion matrix: rows=true class, cols=predicted class
            net_.eval();
            torch::NoGradGuard no_grad;

            for (auto& batch : *test_loader_) {
                auto logits = net_.forward(batch.data);
                auto predicted = logits.argmax(1).to(torch::kCPU);
                auto targets   = batch.target.to(torch::kCPU);

                const int64_t bs = targets.size(0);
                for (int64_t i = 0; i < bs; ++i) {
                    const int64_t t = targets[i].template item<int64_t>();
                    const int64_t p = predicted[i].template item<int64_t>();
                    if (t >= 0 && t < num_classes && p >= 0 && p < num_classes) {
                        (*confusion_matrix_)[static_cast<size_t>(t)][static_cast<size_t>(p)] += 1;
                    }
                }
            }
        }
    }

private:
    NetType& net_;
    torch::optim::Optimizer& optimizer_;
    Loader* train_loader_;
    Loader* test_loader_;
    Loader* val_loader_;
    int epochs_;

    // Logging pointers (owned by caller)
    std::vector<double>* train_loss_ = nullptr;
    std::vector<double>* val_loss_ = nullptr;
    std::vector<double>* test_accuracy_vec_ = nullptr;
    std::vector<double>* val_accuracy_vec_ = nullptr;
    double* test_loss_out_ = nullptr;
    double* test_accuracy_out_ = nullptr;
    std::vector<std::vector<int64_t>>* confusion_matrix_ = nullptr;

    // Computes mean loss per sample + accuracy
    void evaluate(Loader& loader, double* mean_loss_out, double* accuracy_out) {
        net_.eval();
        torch::NoGradGuard no_grad;

        double loss_sum = 0.0;      // sum of per-sample losses
        int64_t correct = 0;
        int64_t total = 0;

        for (auto& batch : loader) {
            auto prediction = net_.forward(batch.data);

            // mean over batch
            auto loss = torch::nn::functional::cross_entropy(prediction, batch.target);

            const int64_t bs = batch.target.size(0);
            loss_sum += loss.template item<double>() * static_cast<double>(bs);

            auto predicted = prediction.argmax(1);
            correct += predicted.eq(batch.target).sum().template item<int64_t>();
            total += bs;
        }

        if (mean_loss_out) {
            *mean_loss_out = (total > 0) ? (loss_sum / static_cast<double>(total)) : 0.0;
        }
        if (accuracy_out) {
            *accuracy_out = (total > 0) ? (static_cast<double>(correct) / static_cast<double>(total)) : 0.0;
        }
    }
};

} // namespace IRL

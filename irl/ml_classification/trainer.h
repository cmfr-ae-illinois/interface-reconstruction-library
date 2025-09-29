#pragma once

#include <torch/torch.h>
#include <vector>
#include <iostream>

namespace IRL {

template<typename NetType, typename Loader>
class Trainer {
public:
    Trainer(NetType& net,
            torch::optim::Optimizer& optimizer,
            Loader* train_loader,
            Loader* test_loader,
            Loader* val_loader,
            int epochs)
        : net_(net),
          optimizer_(optimizer),
          train_loader_(train_loader),
          test_loader_(test_loader),
          val_loader_(val_loader),
          epochs_(epochs) {}

    void train() {
        std::vector<double> lossVector;
        for (int epoch = 1; epoch <= epochs_; ++epoch) {
            net_.train();
            double total_loss = 0.0;
            int correct_predictions = 0;
            int total_samples = 0;

            for (auto& batch : *train_loader_) {
                optimizer_.zero_grad();

                torch::Tensor prediction = net_.forward(batch.data);
                torch::Tensor loss = torch::nn::functional::cross_entropy(prediction, batch.target);
                loss.backward();
                optimizer_.step();

                total_loss += loss.item<double>();
                auto predicted_classes = prediction.argmax(1);
                auto correct = predicted_classes.eq(batch.target);
                correct_predictions += correct.sum().template item<int>();
                total_samples += batch.target.size(0);
            }

            lossVector.push_back(total_loss);

            std::cout << "Epoch [" << epoch << "/" << epochs_ << "] , Accuracy: "
                      << static_cast<double>(correct_predictions) / total_samples << std::endl;

            double test_accuracy = evaluate(*test_loader_);
            std::cout << "Test Accuracy: " << test_accuracy << std::endl;
        }

        double val_accuracy = evaluate(*val_loader_);
        std::cout << "Final Validation Accuracy: " << val_accuracy << std::endl;
    }

private:
    NetType& net_;
    torch::optim::Optimizer& optimizer_;
    Loader* train_loader_;
    Loader* test_loader_;
    Loader* val_loader_;
    int epochs_;

    double evaluate(Loader& loader) {
        net_.eval();
        int correct = 0;
        int total = 0;

        for (auto& batch : loader) {
            auto prediction = net_.forward(batch.data);
            auto predicted = prediction.argmax(1);
            correct += predicted.eq(batch.target).sum().template item<int>();
            total += batch.target.size(0);
        }

        return static_cast<double>(correct) / total;
    }
};

} // namespace IRL

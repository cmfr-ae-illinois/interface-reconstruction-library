#pragma once
#include <torch/torch.h>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace IRL {

class MyDataset : public torch::data::datasets::Dataset<MyDataset> {
private:
    const std::vector<std::vector<float>>* statesV_;
    const std::vector<int>* labelsV_;
    size_t start_;
    size_t end_;
    size_t feature_dim_;

public:
    // Full dataset
    explicit MyDataset(const std::vector<std::vector<float>>* statesV,
                       const std::vector<int>* labelsV)
        : MyDataset(statesV, labelsV, 0, statesV ? statesV->size() : 0) {}

    // Range view: [start, end)
    MyDataset(const std::vector<std::vector<float>>* statesV,
              const std::vector<int>* labelsV,
              size_t start, size_t end)
        : statesV_(statesV),
          labelsV_(labelsV),
          start_(start),
          end_(end),
          feature_dim_(0)
    {
        if (!statesV_ || !labelsV_) {
            throw std::runtime_error("MyDataset: statesV/labelsV is null");
        }
        if (statesV_->size() != labelsV_->size()) {
            throw std::runtime_error("MyDataset: states and labels size mismatch");
        }
        if (start_ > end_ || end_ > statesV_->size()) {
            throw std::runtime_error("MyDataset: invalid range [start,end)");
        }
        if (end_ > start_) {
            feature_dim_ = (*statesV_)[start_].size();
            if (feature_dim_ == 0) {
                throw std::runtime_error("MyDataset: feature dimension is zero");
            }
            // Optional: validate all samples in range have same size (debug safety)
            for (size_t i = start_; i < end_; ++i) {
                if ((*statesV_)[i].size() != feature_dim_) {
                    throw std::runtime_error("MyDataset: inconsistent feature vector size");
                }
            }
        }
    }

    // Return one sample
    torch::data::Example<> get(size_t index) override {
        const size_t i = start_ + index;

        const auto& x = (*statesV_)[i];
        const int y = (*labelsV_)[i];

        // Defensive check: targets for cross_entropy must be in [0, num_classes-1]
        // (You can remove this once you're confident the pipeline is correct.)
        if (y < 0) {
            throw std::runtime_error("MyDataset: negative label encountered: " + std::to_string(y));
        }

        // from_blob would reference vector memory; clone() makes the tensor own its memory (safe for DataLoader workers)
        auto data = torch::from_blob(
                        (void*)x.data(),
                        {(long)x.size()},
                        torch::TensorOptions().dtype(torch::kFloat32)
                    ).clone();

        auto target = torch::tensor((int64_t)y, torch::TensorOptions().dtype(torch::kInt64));

        return {data, target};
    }

    // Number of samples in this view
    torch::optional<size_t> size() const override {
        return (end_ >= start_) ? (end_ - start_) : 0;
    }

    // Optional helper
    size_t feature_dim() const { return feature_dim_; }
};

} // namespace IRL


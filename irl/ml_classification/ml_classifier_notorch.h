#pragma once

#include "classifier.h"

#include <vector>

namespace IRL {

class MLClassifierNoTorch : public Classifier {
public:
    MLClassifierNoTorch();

    int classify(const std::vector<float>& flattened_state,
                 std::vector<float>* out_probs = nullptr) override;

private:
    void forwardLogits(const std::vector<float>& flattened_state,
                       std::vector<double>& logits) const;

    static std::vector<float> stableSoftmax(const std::vector<double>& logits);

    static int argmax(const std::vector<double>& values);
};

}  // namespace IRL

#include "ml_classifier_notorch.tpp"
#pragma once

#include "ml_classifier_weights_and_biases.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace IRL {

inline MLClassifierNoTorch::MLClassifierNoTorch(int stencil_size)
    : Classifier(stencil_size) {}

inline void MLClassifierNoTorch::forwardLogits(
    const std::vector<float>& flattened_state,
    std::vector<double>& logits) const {
    if (static_cast<int>(flattened_state.size()) != mlclassifier::input_size) {
        throw std::runtime_error(
            "MLClassifierNoTorch::forwardLogits(): flattened_state has wrong size.");
    }

    for (float v : flattened_state) {
        if (!std::isfinite(v)) {
            throw std::runtime_error(
                "MLClassifierNoTorch::forwardLogits(): flattened_state contains NaN or Inf.");
        }
    }

    std::array<double, mlclassifier::hidden_size1> h1{};
    std::array<double, mlclassifier::hidden_size2> h2{};
    std::array<double, mlclassifier::hidden_size3> h3{};

    for (int j = 0; j < mlclassifier::hidden_size1; ++j) {
        h1[j] = mlclassifier::fc1_bias[j];

        for (int i = 0; i < mlclassifier::input_size; ++i) {
            h1[j] += static_cast<double>(flattened_state[i]) *
                     mlclassifier::fc1_weight[j][i];
        }

        if (h1[j] < 0.0) {
            h1[j] = 0.0;
        }
    }

    for (int j = 0; j < mlclassifier::hidden_size2; ++j) {
        h2[j] = mlclassifier::fc2_bias[j];

        for (int i = 0; i < mlclassifier::hidden_size1; ++i) {
            h2[j] += h1[i] * mlclassifier::fc2_weight[j][i];
        }

        if (h2[j] < 0.0) {
            h2[j] = 0.0;
        }
    }

    for (int j = 0; j < mlclassifier::hidden_size3; ++j) {
        h3[j] = mlclassifier::fc3_bias[j];

        for (int i = 0; i < mlclassifier::hidden_size2; ++i) {
            h3[j] += h2[i] * mlclassifier::fc3_weight[j][i];
        }

        if (h3[j] < 0.0) {
            h3[j] = 0.0;
        }
    }

    logits.assign(mlclassifier::output_size, 0.0);

    for (int j = 0; j < mlclassifier::output_size; ++j) {
        logits[j] = mlclassifier::fc4_bias[j];

        for (int i = 0; i < mlclassifier::hidden_size3; ++i) {
            logits[j] += h3[i] * mlclassifier::fc4_weight[j][i];
        }
    }
}

inline int MLClassifierNoTorch::classify(
    const std::vector<float>& flattened_state,
    std::vector<float>* out_probs) {
    std::vector<double> logits;
    forwardLogits(flattened_state, logits);

    const int predicted_class = argmax(logits);

    if (out_probs) {
        *out_probs = stableSoftmax(logits);
    }

    return predicted_class;
}

inline std::vector<float> MLClassifierNoTorch::stableSoftmax(
    const std::vector<double>& logits) {
    if (logits.empty()) {
        throw std::runtime_error(
            "MLClassifierNoTorch::stableSoftmax(): empty logits.");
    }

    const double max_logit =
        *std::max_element(logits.begin(), logits.end());

    std::vector<float> probs(logits.size(), 0.0f);
    double sum = 0.0;

    for (std::size_t i = 0; i < logits.size(); ++i) {
        const double e = std::exp(logits[i] - max_logit);
        probs[i] = static_cast<float>(e);
        sum += e;
    }

    if (!(sum > 0.0) || !std::isfinite(sum)) {
        throw std::runtime_error(
            "MLClassifierNoTorch::stableSoftmax(): invalid softmax denominator.");
    }

    for (float& p : probs) {
        p = static_cast<float>(static_cast<double>(p) / sum);
    }

    return probs;
}

inline int MLClassifierNoTorch::argmax(const std::vector<double>& values) {
    if (values.empty()) {
        throw std::runtime_error(
            "MLClassifierNoTorch::argmax(): empty vector.");
    }

    return static_cast<int>(
        std::distance(values.begin(),
                      std::max_element(values.begin(), values.end())));
}

}  // namespace IRL
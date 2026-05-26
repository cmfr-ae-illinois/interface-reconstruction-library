#pragma once

#include "classifier.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace IRL {

class MLClassifierNoTorch : public Classifier {
private:
    int stencil_size = 3;
    int input_size = 27;
    int hidden_size1 = 128;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 4;

    bool weights_loaded = false;

    // Row-major matrices:
    // fc1_weight[out_index * input_size + in_index]
    std::vector<float> fc1_weight;
    std::vector<float> fc1_bias;

    std::vector<float> fc2_weight;
    std::vector<float> fc2_bias;

    std::vector<float> fc3_weight;
    std::vector<float> fc3_bias;

    std::vector<float> fc4_weight;
    std::vector<float> fc4_bias;

public:
    MLClassifierNoTorch(int stencil = 3,
                 int input = 27,
                 int h1 = 128,
                 int h2 = 64,
                 int h3 = 32,
                 int out = 4)
        : Classifier(stencil),
          stencil_size(stencil),
          input_size(input),
          hidden_size1(h1),
          hidden_size2(h2),
          hidden_size3(h3),
          output_size(out) {
        allocateWeights();
    }

    MLClassifierNoTorch(const std::string& weights_file,
                 int stencil = 3,
                 int input = 27,
                 int h1 = 128,
                 int h2 = 64,
                 int h3 = 32,
                 int out = 4)
        : MLClassifierNoTorch(stencil, input, h1, h2, h3, out) {
        loadWeights(weights_file);
    }

    void loadWeights(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        if (!in) {
            throw std::runtime_error("MLClassifier::loadWeights(): failed to open " + filename);
        }

        constexpr char expected_magic[8] = {'I', 'R', 'L', 'M', 'L', 'P', '1', '\0'};

        char magic[8] = {};
        readExact(in, magic, sizeof(magic), "magic");

        if (std::memcmp(magic, expected_magic, sizeof(magic)) != 0) {
            throw std::runtime_error(
                "MLClassifier::loadWeights(): wrong file format. "
                "Expected IRLMLP1 binary weights file.");
        }

        std::uint32_t version = 0;
        readPod(in, version, "version");

        if (version != 1) {
            throw std::runtime_error("MLClassifier::loadWeights(): unsupported weights file version.");
        }

        std::int32_t file_input = 0;
        std::int32_t file_h1 = 0;
        std::int32_t file_h2 = 0;
        std::int32_t file_h3 = 0;
        std::int32_t file_out = 0;

        readPod(in, file_input, "input_size");
        readPod(in, file_h1, "hidden_size1");
        readPod(in, file_h2, "hidden_size2");
        readPod(in, file_h3, "hidden_size3");
        readPod(in, file_out, "output_size");

        if (file_input != input_size ||
            file_h1 != hidden_size1 ||
            file_h2 != hidden_size2 ||
            file_h3 != hidden_size3 ||
            file_out != output_size) {
            throw std::runtime_error(
                "MLClassifier::loadWeights(): network dimensions in file do not match "
                "the MLClassifier constructor arguments.");
        }

        readVector(in, fc1_weight, checkedMatrixSize(hidden_size1, input_size), "fc1_weight");
        readVector(in, fc1_bias, checkedSize(hidden_size1), "fc1_bias");

        readVector(in, fc2_weight, checkedMatrixSize(hidden_size2, hidden_size1), "fc2_weight");
        readVector(in, fc2_bias, checkedSize(hidden_size2), "fc2_bias");

        readVector(in, fc3_weight, checkedMatrixSize(hidden_size3, hidden_size2), "fc3_weight");
        readVector(in, fc3_bias, checkedSize(hidden_size3), "fc3_bias");

        readVector(in, fc4_weight, checkedMatrixSize(output_size, hidden_size3), "fc4_weight");
        readVector(in, fc4_bias, checkedSize(output_size), "fc4_bias");

        weights_loaded = true;
    }

    std::vector<float> forwardLogits(const std::vector<float>& flattened_state) const {
        if (!weights_loaded) {
            throw std::runtime_error("MLClassifier::forwardLogits(): weights have not been loaded.");
        }

        if (static_cast<int>(flattened_state.size()) != input_size) {
            throw std::runtime_error(
                "MLClassifier::forwardLogits(): flattened_state has wrong size.");
        }

        for (float v : flattened_state) {
            if (!std::isfinite(v)) {
                throw std::runtime_error(
                    "MLClassifier::forwardLogits(): flattened_state contains NaN or Inf.");
            }
        }

        std::vector<float> h1;
        std::vector<float> h2;
        std::vector<float> h3;
        std::vector<float> logits;

        linear(flattened_state, fc1_weight, fc1_bias, input_size, hidden_size1, h1, true);
        linear(h1, fc2_weight, fc2_bias, hidden_size1, hidden_size2, h2, true);
        linear(h2, fc3_weight, fc3_bias, hidden_size2, hidden_size3, h3, true);
        linear(h3, fc4_weight, fc4_bias, hidden_size3, output_size, logits, false);

        return logits;
    }

    int classify(const std::vector<float>& flattened_state,
                 std::vector<float>* out_probs = nullptr) override {
        std::vector<float> logits = forwardLogits(flattened_state);

        const int predicted_class = argmax(logits);

        if (out_probs) {
            *out_probs = stableSoftmax(logits);
        }

        return predicted_class;
    }

private:
    void allocateWeights() {
        fc1_weight.resize(checkedMatrixSize(hidden_size1, input_size));
        fc1_bias.resize(checkedSize(hidden_size1));

        fc2_weight.resize(checkedMatrixSize(hidden_size2, hidden_size1));
        fc2_bias.resize(checkedSize(hidden_size2));

        fc3_weight.resize(checkedMatrixSize(hidden_size3, hidden_size2));
        fc3_bias.resize(checkedSize(hidden_size3));

        fc4_weight.resize(checkedMatrixSize(output_size, hidden_size3));
        fc4_bias.resize(checkedSize(output_size));
    }

    static std::size_t checkedSize(int n) {
        if (n <= 0) {
            throw std::runtime_error("MLClassifier: invalid non-positive layer size.");
        }
        return static_cast<std::size_t>(n);
    }

    static std::size_t checkedMatrixSize(int rows, int cols) {
        if (rows <= 0 || cols <= 0) {
            throw std::runtime_error("MLClassifier: invalid non-positive matrix size.");
        }

        return static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    }

    static void readExact(std::ifstream& in,
                          char* dst,
                          std::size_t nbytes,
                          const std::string& what) {
        in.read(dst, static_cast<std::streamsize>(nbytes));

        if (!in) {
            throw std::runtime_error("MLClassifier::loadWeights(): failed while reading " + what);
        }
    }

    template <typename T>
    static void readPod(std::ifstream& in, T& value, const std::string& what) {
        readExact(in, reinterpret_cast<char*>(&value), sizeof(T), what);
    }

    static void readVector(std::ifstream& in,
                           std::vector<float>& dst,
                           std::size_t count,
                           const std::string& what) {
        dst.resize(count);
        readExact(in,
                  reinterpret_cast<char*>(dst.data()),
                  count * sizeof(float),
                  what);
    }

    static void linear(const std::vector<float>& x,
                       const std::vector<float>& weight,
                       const std::vector<float>& bias,
                       int in_features,
                       int out_features,
                       std::vector<float>& y,
                       bool apply_relu) {
        y.assign(static_cast<std::size_t>(out_features), 0.0f);

        for (int out = 0; out < out_features; ++out) {
            double acc = static_cast<double>(bias[static_cast<std::size_t>(out)]);

            const std::size_t row_offset =
                static_cast<std::size_t>(out) * static_cast<std::size_t>(in_features);

            for (int in = 0; in < in_features; ++in) {
                acc += static_cast<double>(weight[row_offset + static_cast<std::size_t>(in)]) *
                       static_cast<double>(x[static_cast<std::size_t>(in)]);
            }

            if (apply_relu && acc < 0.0) {
                acc = 0.0;
            }

            y[static_cast<std::size_t>(out)] = static_cast<float>(acc);
        }
    }

    static std::vector<float> stableSoftmax(const std::vector<float>& logits) {
        if (logits.empty()) {
            throw std::runtime_error("MLClassifier::stableSoftmax(): empty logits.");
        }

        const float max_logit = *std::max_element(logits.begin(), logits.end());

        std::vector<float> probs(logits.size(), 0.0f);
        double sum = 0.0;

        for (std::size_t i = 0; i < logits.size(); ++i) {
            const double e = std::exp(static_cast<double>(logits[i] - max_logit));
            probs[i] = static_cast<float>(e);
            sum += e;
        }

        if (!(sum > 0.0) || !std::isfinite(sum)) {
            throw std::runtime_error("MLClassifier::stableSoftmax(): invalid softmax denominator.");
        }

        for (float& p : probs) {
            p = static_cast<float>(static_cast<double>(p) / sum);
        }

        return probs;
    }

    static int argmax(const std::vector<float>& values) {
        if (values.empty()) {
            throw std::runtime_error("MLClassifier::argmax(): empty vector.");
        }

        return static_cast<int>(
            std::distance(values.begin(),
                          std::max_element(values.begin(), values.end())));
    }
};

} // namespace IRL
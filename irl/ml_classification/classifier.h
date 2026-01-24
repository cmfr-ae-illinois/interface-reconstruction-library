#pragma once
#include <vector>

namespace IRL {

// Abstract base classifier
class Classifier {
protected:
    int stencil_size;

public:
    explicit Classifier(int stencil) : stencil_size(stencil) {}
    virtual ~Classifier() = default;

    int getStencilSize() const { return stencil_size; }

    // Every classifier must implement classify
    virtual int classify(const std::vector<float>& flattened_state, std::vector<float>* out_probs = nullptr) = 0;
};

} // namespace IRL
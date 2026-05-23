#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

#include <pybind11/pybind11.h>

namespace py = pybind11;

int longest_common_substring_len(const std::string& a, const std::string& b) {
    if (a.empty() || b.empty()) {
        return 0;
    }

    std::vector<int> prev(b.size() + 1, 0);
    int best = 0;
    for (char ca : a) {
        std::vector<int> curr(b.size() + 1, 0);
        for (size_t j = 0; j < b.size(); ++j) {
            if (ca == b[j]) {
                curr[j + 1] = prev[j] + 1;
                if (curr[j + 1] > best) {
                    best = curr[j + 1];
                }
            }
        }
        prev.swap(curr);
    }
    return best;
}

std::tuple<std::string, std::string, int> needleman_wunsch(
    const std::string& ref,
    const std::string& seq,
    int match = 2,
    int mismatch = -1,
    int gap = -2
) {
    const size_t n = ref.size();
    const size_t m = seq.size();

    std::vector<int> score((n + 1) * (m + 1), 0);
    std::vector<unsigned char> trace((n + 1) * (m + 1), 0);

    auto idx = [m](size_t i, size_t j) {
        return i * (m + 1) + j;
    };

    for (size_t i = 1; i <= n; ++i) {
        score[idx(i, 0)] = score[idx(i - 1, 0)] + gap;
        trace[idx(i, 0)] = 1;
    }
    for (size_t j = 1; j <= m; ++j) {
        score[idx(0, j)] = score[idx(0, j - 1)] + gap;
        trace[idx(0, j)] = 2;
    }

    for (size_t i = 1; i <= n; ++i) {
        for (size_t j = 1; j <= m; ++j) {
            int diag = score[idx(i - 1, j - 1)] + (ref[i - 1] == seq[j - 1] ? match : mismatch);
            int up = score[idx(i - 1, j)] + gap;
            int left = score[idx(i, j - 1)] + gap;

            int best = diag;
            unsigned char move = 0;
            if (up > best) {
                best = up;
                move = 1;
            }
            if (left > best) {
                best = left;
                move = 2;
            }

            score[idx(i, j)] = best;
            trace[idx(i, j)] = move;
        }
    }

    size_t i = n;
    size_t j = m;
    std::string aligned_ref;
    std::string aligned_seq;
    aligned_ref.reserve(n + m);
    aligned_seq.reserve(n + m);

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && trace[idx(i, j)] == 0) {
            aligned_ref.push_back(ref[i - 1]);
            aligned_seq.push_back(seq[j - 1]);
            --i;
            --j;
        } else if (i > 0 && (j == 0 || trace[idx(i, j)] == 1)) {
            aligned_ref.push_back(ref[i - 1]);
            aligned_seq.push_back('-');
            --i;
        } else {
            aligned_ref.push_back('-');
            aligned_seq.push_back(seq[j - 1]);
            --j;
        }
    }

    std::reverse(aligned_ref.begin(), aligned_ref.end());
    std::reverse(aligned_seq.begin(), aligned_seq.end());

    return std::make_tuple(aligned_ref, aligned_seq, score[idx(n, m)]);
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "Fast C++ alignment helpers for SCAN_SPLASH";
    m.def(
        "longest_common_substring_len",
        &longest_common_substring_len,
        py::arg("a"),
        py::arg("b")
    );
    m.def(
        "needleman_wunsch",
        &needleman_wunsch,
        py::arg("ref"),
        py::arg("seq"),
        py::arg("match") = 2,
        py::arg("mismatch") = -1,
        py::arg("gap") = -2
    );
}

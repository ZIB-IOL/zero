// Independent complete toroidal label audit. Standard library only.
// c++ -O3 -std=c++17 independent_labels.cpp -o independent_labels
// independent_labels NX NY NZ labels.bin exact_mask.csv
// Rows use std::bitset<256>, including arbitrary third-axis lengths <=256.
#include <array>
#include <bitset>
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
        if (argc != 6) throw std::runtime_error("usage: independent_labels NX NY NZ labels.bin exact_mask.csv");
        auto started = std::chrono::steady_clock::now();
        std::array<int, 3> shape{std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3])};
        if (shape[0] <= 0 || shape[1] <= 0 || shape[2] <= 0 || shape[2] > 256)
            throw std::runtime_error("invalid shape (third axis must be <=256)");
        size_t rows = size_t(shape[0]) * shape[1], total = rows * shape[2];
        using Bits = std::bitset<256>;
        std::vector<std::array<Bits,14>> color_rows(rows);
        std::vector<unsigned> present(rows, 0);
        std::array<size_t,15> counts{};
        std::ifstream labels(argv[4], std::ios::binary);
        if (!labels) throw std::runtime_error("cannot open labels");
        for (size_t index = 0; index < total; ++index) {
            int color = labels.get();
            if (color < 0 || color > 14) throw std::runtime_error("invalid label or truncated labels");
            ++counts[color];
            if (color == 14) continue;
            size_t row = index / shape[2], position = index % shape[2];
            color_rows[row][color].set(position);
            present[row] |= 1u << color;
        }
        if (labels.get() != EOF) throw std::runtime_error("extra bytes in labels");
        for (int c = 0; c < 14; ++c)
            if (!counts[c]) throw std::runtime_error("not exactly fourteen real colors");
        std::ifstream mask(argv[5]);
        if (!mask) throw std::runtime_error("cannot open exact mask");
        std::set<std::array<int,3>> offsets;
        std::string line;
        while (std::getline(mask,line)) {
            if (line.empty()) throw std::runtime_error("empty mask row");
            std::istringstream parser(line);
            std::array<int,3> d{};
            char comma1, comma2;
            if (!(parser >> d[0] >> comma1 >> d[1] >> comma2 >> d[2]) || comma1 != ',' || comma2 != ',')
                throw std::runtime_error("invalid mask row");
            parser >> std::ws;
            if (!parser.eof()) throw std::runtime_error("extra mask columns");
            for (int axis = 0; axis < 3; ++axis)
                d[axis] = (d[axis] % shape[axis] + shape[axis]) % shape[axis];
            offsets.insert(d);
        }
        if (offsets.empty()) throw std::runtime_error("empty mask");
        size_t violations = 0, checked = 0;
        for (const auto& d : offsets) {
            for (int i = 0; i < shape[0]; ++i) {
                int next_i = (i + d[0]) % shape[0];
                for (int j = 0; j < shape[1]; ++j) {
                    size_t first = size_t(i)*shape[1] + j;
                    size_t second = size_t(next_i)*shape[1] + (j+d[1])%shape[1];
                    unsigned shared = present[first] & present[second];
                    while (shared) {
                        unsigned color = __builtin_ctz(shared);
                        shared &= shared - 1;
                        const Bits& source = color_rows[second][color];
                        // Right rotation within NZ bits. Bits above NZ are zero
                        // in first, so the final intersection discards spillover.
                        Bits shifted = source >> d[2];
                        if (d[2]) shifted |= source << (shape[2] - d[2]);
                        Bits hits = shifted & color_rows[first][color];
                        if (hits.any()) {
                            violations += hits.count();
                            if (violations < 10)
                                std::cerr << "conflict row=" << first << " color=" << color
                                          << " shift=" << d[0] << ',' << d[1] << ',' << d[2] << '\n';
                        }
                    }
                }
            }
            if (++checked % 10000 == 0)
                std::cerr << "independent offsets=" << checked << " violations=" << violations << '\n';
        }
        double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        std::cout << "{\"verified\":" << (violations ? "false" : "true")
                  << ",\"audited_offsets\":" << checked << ",\"monochromatic_conflict_pairs\":" << violations
                  << ",\"total_cells\":" << total << ",\"removed_count\":" << counts[14]
                  << ",\"color_counts\":[";
        for (int c = 0; c < 15; ++c) std::cout << (c ? "," : "") << counts[c];
        std::cout << "],\"seconds\":" << seconds << "}\n";
        return violations ? 2 : 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

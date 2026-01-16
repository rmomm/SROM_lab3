#pragma once
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>

class GF2_359 {
public:
    static constexpr int m = 359;
    static constexpr int w = 6; 
    std::array<uint64_t, w> data;

    GF2_359();
    explicit GF2_359(uint64_t v);

    static GF2_359 zero();
    static GF2_359 one();

    GF2_359 operator+(const GF2_359& other) const;
    GF2_359& operator+=(const GF2_359& other);

    GF2_359 operator*(const GF2_359& other) const;
    GF2_359& operator*=(const GF2_359& other);
    GF2_359 square() const;
    GF2_359 pow(uint64_t e) const;
    GF2_359 inverse() const;

    bool operator==(const GF2_359& other) const;
    bool operator!=(const GF2_359& other) const;

    bool get_bit(int i) const;
    void set_bit(int i, bool v);

    int degree() const;
    GF2_359 shift_left(int k) const;
    void normalize();

    int trace() const;

    const std::array<uint64_t, w>& get_data() const;
    std::string toHex() const;
    static GF2_359 fromHex(const std::string& hex);

private:
    GF2_359 mul_reduce(const GF2_359& a, const GF2_359& b) const;
    void reduce(uint64_t temp[12]) const;
};

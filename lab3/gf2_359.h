#pragma once
#include <vector>
#include <string>
#include <cstdint>

using namespace std;

class GF2_359 {
public:
    GF2_359();
    explicit GF2_359(const vector<uint64_t>& data);

    static GF2_359 fromHex(const string& hex);
    string toHex() const;

    static GF2_359 zero();
    static GF2_359 one();

    GF2_359 operator+(const GF2_359& other) const;
    GF2_359 operator*(const GF2_359& other) const;

    GF2_359 square() const;
    GF2_359 inverse() const;
    GF2_359 pow(const GF2_359& exponent) const;
    int trace() const;

    bool operator==(const GF2_359& other) const;
    bool operator!=(const GF2_359& other) const;

private:
    vector<uint64_t> bits;

    static int degree(const vector<uint64_t>& a);
    static void trim(vector<uint64_t>& a);
    static bool isOne(const vector<uint64_t>& a);

    static vector<uint64_t> xorPoly(
        const vector<uint64_t>& a,
        const vector<uint64_t>& b
    );

    static vector<uint64_t> shiftLeft(
        const vector<uint64_t>& a,
        int shift
    );

    static void divMod(
        const vector<uint64_t>& a,
        const vector<uint64_t>& b,
        vector<uint64_t>& q,
        vector<uint64_t>& r
    );

    static vector<uint64_t> mulPoly(
        const vector<uint64_t>& a,
        const vector<uint64_t>& b
    );

    static vector<uint64_t> reduce(
        const vector<uint64_t>& a
    );
};


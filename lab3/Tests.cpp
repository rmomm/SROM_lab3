#include "gf2_359.h"
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace chrono;

void TestGF2m() {
    GF2_359 A = GF2_359::fromHex("9f3c8a71d4b29e6c8b5fa31d27c4e9a6");
    GF2_359 B = GF2_359::fromHex("7d1a9e8f63b0c4d52a79f8e6b3d2c1a0");
    GF2_359 C = GF2_359::fromHex("abcdef1234567890fedcba9876543210");

    cout << "A = " << A.toHex() << endl;
    cout << "B = " << B.toHex() << endl;
    cout << "C = " << C.toHex() << endl;

    GF2_359 R;

    R = A + B;
    cout << "A + B = " << R.toHex() << endl;

    R = A * B;
    cout << "A * B = " << R.toHex() << endl;

    R = A.square();
    cout << "A^2 = " << R.toHex() << endl;

    R = A.pow(5);
    cout << "A^5 = " << R.toHex() << endl;

    R = A.inverse();
    cout << "A^-1 = " << R.toHex() << endl;

    R = (A + B) * C;
    GF2_359 rhs = A * C + B * C;
    cout << "(A+B)*C = " << R.toHex() << endl;
    cout << "A*C + B*C = " << rhs.toHex() << endl;
    cout << "Distributive law holds? " << (R == rhs ? "YES" : "NO") << endl;
}

void BenchmarkGF2m() {
    const int iter = 1000;

    GF2_359 A = GF2_359::fromHex("9f3c8a71d4b29e6c8b5fa31d27c4e9a6");
    GF2_359 B = GF2_359::fromHex("7d1a9e8f63b0c4d52a79f8e6b3d2c1a0");
    GF2_359 C = GF2_359::fromHex("abcdef1234567890fedcba9876543210");
    GF2_359 R;

    using clock = high_resolution_clock;

    auto bench = [&](const string& name, auto func) {
        auto t1 = clock::now();
        for (int i = 0; i < iter; ++i) func();
        auto t2 = clock::now();
        long long avg = chrono::duration_cast<nanoseconds>(t2 - t1).count() / iter;
        cout << setw(15) << left << name << avg << " ns" << endl;
        };

    cout << "\nAverage time per operation (GF2_359):\n";
    cout << "-------------------------------------\n";

    bench("add", [&]() { R = A + B; });
    bench("mul", [&]() { R = A * B; });
    bench("square", [&]() { R = A.square(); });
    bench("pow 5", [&]() { R = A.pow(5); });
    bench("inverse", [&]() { R = A.inverse(); });

}

int main() {
    TestGF2m();
    BenchmarkGF2m();
    return 0;
}

#include <cassert>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <array>
#include <vector>
#include <string>
#include <complex>
#include <cmath>

const long double PI = acos(-1);
using complex = std::complex<long double>;
using std::cout;
using std::cerr;
using std::endl;

struct BigInteger {
private:
    static const long long BASE = 10000;
    static const int DIGIT_SIZE = 4;

    bool negative = false; // Если равно 0, тоже false
    std::vector<long long> digit; // "Цифры" хранятся - digit[i] * BASE**i, без лидирующих нулей

    friend class Rational;

    static void fft(std::vector<complex>& a, bool inverse) {
        size_t n = a.size();
        for (size_t i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                std::swap(a[i], a[j]);
        }
        for (size_t len = 2; len <= n; len *= 2) {
            long double angle = 2 * PI / static_cast<long double> (len) * (inverse ? -1 : 1);
            complex primitive_root(cosl(angle), sinl(angle));
            for (size_t i = 0; i < n; i += len) {
                complex w(1);
                complex* ptrL = &a[i];
                complex* ptrR = ptrL + len / 2;
                complex* ptrEnd = ptrL + len;
                for (; ptrR != ptrEnd; ++ptrL, ++ptrR) {
                    complex u = *ptrL;
                    complex v = *ptrR * w;
                    *ptrL = u + v;
                    *ptrR = u - v;
                    w *= primitive_root;
                }
            }
        }
        if (inverse) {
            for (complex& x: a) {
                x /= n;
            }
        }
    }

    // Прибавляет к abs(*this) abs(x) * sign, причем abs(*this) гарантированно больше abs(x)
    void _add(const BigInteger& x, bool sign) {
        long long carry = 0;
        for (size_t i = 0; i < std::min(amountDigits(), x.amountDigits()) || carry; ++i) {
            if (i >= amountDigits()) {
                digit.push_back(0);
            }
            digit[i] += (sign ? -1 : 1) * carry;
            if (i < x.amountDigits()) {
                digit[i] += (sign ? -1 : 1) * x.digit[i];
            }
            carry = 0;
            if (digit[i] >= BASE) {
                carry = 1;
                digit[i] -= BASE;
            }
            if (digit[i] < 0) {
                carry = 1;
                digit[i] += BASE;
            }
        }
        normalize();
    }

    [[nodiscard]] BigInteger abs() const;

    // Удаляет лидирующие нули
    void normalize() {
        while (amountDigits() > 1 && digit.back() == 0)
            digit.pop_back();
        if (amountDigits() == 1 && digit.back() == 0)
            negative = false;
    }

public:
    friend std::istream& operator>>(std::istream& in, BigInteger& x);
    BigInteger& operator*=(int y);
    BigInteger& operator+=(const BigInteger& x);
    BigInteger& operator-=(const BigInteger& x);
    BigInteger& operator++();
    BigInteger& operator--();
    BigInteger& operator*=(const BigInteger& y);
    BigInteger operator*(const BigInteger& y) const;
    BigInteger operator/(const BigInteger& y) const;
    BigInteger operator%(const BigInteger& y) const;
    BigInteger& operator/=(BigInteger y);
    BigInteger& operator/=(int y);
    BigInteger& operator%=(const BigInteger& y);
    const BigInteger operator++(int);
    const BigInteger operator--(int);
    void swap(BigInteger& oth);

    template<typename T>
    BigInteger(T x) {
        if (x < 0) {
            x *= -1;
            negative = true;
        }
        do {
            digit.push_back(x % BASE);
            x /= BASE;
        } while (x != 0);
    }

    BigInteger(const BigInteger& oth) = default;

    BigInteger() : digit({0}) {}

    // Возвращает количество "цифр"
    [[nodiscard]] size_t amountDigits() const {
        return digit.size();
    }

    // true, если положительно
    [[nodiscard]] bool sign() const {
        return !negative;
    }

    ~BigInteger() {
        digit.clear();
        digit.shrink_to_fit();
    }

    explicit operator bool() const {
        return amountDigits() > 1 || digit[0] != 0;
    }

    bool operator<(const BigInteger& oth) const {
        if (sign() != oth.sign())
            return !sign();
        if (amountDigits() != oth.amountDigits()) {
            return (sign()
                ^ (amountDigits() > oth.amountDigits())); // Одно из двух: либо мы положительны, либо у нас больше цифр
        }
        for (int i = static_cast <int> (amountDigits()) - 1; i >= 0; --i) {
            if (digit[i] != oth.digit[i]) {
                return sign() ^ (digit[i] > oth.digit[i]); // Одно из двух: либо мы положительны, либо наша цифра больше
            }
        }
        return false;
    }

    bool operator==(const BigInteger& oth) const {
        return sign() == oth.sign() && digit == oth.digit;
    }

    bool operator!=(const BigInteger& oth) const {
        return !(*this == oth);
    }

    bool operator>(const BigInteger& oth) const {
        return oth < *this;
    }

    bool operator<=(const BigInteger& oth) const {
        return !(oth < *this);
    }

    bool operator>=(const BigInteger& oth) const {
        return !(*this < oth);
    }

    BigInteger& operator=(const BigInteger& oth)& {
        if (this == &oth)
            return *this;
        BigInteger y = oth;
        this->swap(y); // Просто swap(y) выглядит очень странно :)
        return *this;
    }

    BigInteger operator-() const {
        BigInteger x(*this);
        x.negative ^= 1;
        x.normalize(); // Пусть будет
        return x;
    }

    [[nodiscard]] std::string toString() const {
        std::string res;
        if (negative)
            res += '-';
        for (int i = static_cast<int>(amountDigits()) - 1; i >= 0; --i) {
            std::string currentDigit = std::to_string(digit[i]);
            if (i != static_cast<int>(amountDigits()) - 1) {
                res += std::string(DIGIT_SIZE - currentDigit.size(), '0'); // Добавляем опущенные нули
            }
            res += currentDigit;
        }
        return res;
    }

};

void BigInteger::swap(BigInteger& oth) {
    digit.swap(oth.digit);
    std::swap(oth.negative, negative);
}

std::ostream& operator<<(std::ostream& out, const BigInteger& x) {
    out << x.toString();
    return out;
}

BigInteger BigInteger::abs() const {
    BigInteger ans(*this);
    ans.negative = false;
    return ans;
}

BigInteger operator "" _bi(unsigned long long x) {
    return BigInteger(x);
}

std::istream& operator>>(std::istream& in, BigInteger& x) {
    std::string s;
    in >> s;
    x.digit.clear();
    int firstDigitPos = 0;
    x.negative = false;
    if (s[0] == '-') {
        x.negative = true;
        firstDigitPos = 1;
    } else if (s[0] == '+') {
        firstDigitPos = 1;
    }
    for (int digitEnd = static_cast<int>(s.size()); digitEnd > firstDigitPos; digitEnd -= BigInteger::DIGIT_SIZE) {
        if (digitEnd - BigInteger::DIGIT_SIZE >= firstDigitPos)
            x.digit.push_back(std::stoll(s.substr(digitEnd - BigInteger::DIGIT_SIZE, BigInteger::DIGIT_SIZE)));
        else
            x.digit.push_back(std::stoll(s.substr(firstDigitPos, digitEnd - firstDigitPos)));
    }
    if (!x.amountDigits()) {
        x.digit.push_back(0);
    }
    x.normalize();
    return in;
}

BigInteger operator-(const BigInteger& x, const BigInteger& y);
BigInteger operator+(const BigInteger& x, const BigInteger& y);

BigInteger& BigInteger::operator-=(const BigInteger& x) {
    return *this = *this + (-x);
}

BigInteger operator-(const BigInteger& x, const BigInteger& y) {
    BigInteger ans = x;
    ans -= y;
    return ans;
}

BigInteger& BigInteger::operator+=(const BigInteger& x) {
    if (this->abs() < x.abs()) {
        return *this = x + *this;
    }
    _add(x, negative ^ x.negative);
    return *this;
}

BigInteger operator+(const BigInteger& x, const BigInteger& y) {
    BigInteger ans = x;
    ans += y;
    return ans;
}

BigInteger& BigInteger::operator++() {
    return (*this += 1_bi);
}

BigInteger& BigInteger::operator--() {
    return (*this -= 1_bi);
}

const BigInteger BigInteger::operator++(int) {
    BigInteger res(*this);
    *this += 1_bi;
    return res;
}

const BigInteger BigInteger::operator--(int) {
    BigInteger res = *this;
    *this -= 1_bi;
    return res;
}

BigInteger& BigInteger::operator*=(const BigInteger& y) {
    negative ^= y.negative;
    std::vector<complex> fx(digit.begin(), digit.end());
    std::vector<complex> fy(y.digit.begin(), y.digit.end());
    size_t n = 1;
    while (n < std::max(fx.size(), fy.size())) {
        n *= 2;
    }
    n *= 2;
    fx.resize(n);
    fy.resize(n);
    fft(fx, false);
    fft(fy, false);
    for (size_t i = 0; i < n; ++i) {
        fx[i] *= fy[i];
    }
    fft(fx, true);
    long long carry = 0;
    for (size_t i = 0; i < n || carry; ++i) {
        if (i >= amountDigits()) {
            digit.push_back(0);
        }
        long long cur = carry;
        if (i < fx.size())
            cur += llroundl(fx[i].real());
        carry = 0;
        if (cur >= BASE) {
            carry = cur / BASE;
            cur %= BASE;
        }
        digit[i] = cur;
    }
    normalize();
    return *this;
}

BigInteger& BigInteger::operator*=(int y) {
    if (llabs(y) >= BASE) {
        return *this *= BigInteger(y);
    }
    if (y < 0) {
        negative ^= 1;
        y = -y;
    }
    int carry = 0;
    for (size_t i = 0; i < digit.size() || carry; ++i) {
        if (i == digit.size())
            digit.push_back(0);
        long long cur = digit[i] * static_cast <long long>(y) + carry;
        carry = static_cast <int> (cur / BASE);
        digit[i] = (cur % BASE);
    }
    normalize();
    return *this;
}

BigInteger BigInteger::operator*(const BigInteger& y) const {
    BigInteger res(*this);
    return (res *= y);
}

BigInteger& BigInteger::operator/=(int y) {
    long long carry = 0;
    for (int i = (int) amountDigits() - 1; i >= 0; --i) {
        long long cur = digit[i] + carry * BASE;
        digit[i] = cur / y;
        carry = cur % y;
    }
    normalize();
    return *this;
}

BigInteger BigInteger::operator/(const BigInteger& y) const {
    BigInteger res(*this);
    return (res /= y);
}

BigInteger& BigInteger::operator/=(BigInteger y) {
    auto cl = clock();
    bool f = false;
    if (digit.size() > 200 || y.digit.size() > 200) {
        f = true;
        std::cerr << "division " << digit.size() << ' ' << y.digit.size() << std::endl;
    }
    bool resultingNegativeness = negative ^ y.negative;
    negative = false;
    y.negative = false;
    if (*this < y) {
        *this = 0;
        return *this;
    }
    std::vector<int> digitsBigEndian(amountDigits() * DIGIT_SIZE);
    int* ptr = &digitsBigEndian[0];
    for (long long dig: digit) {
        for (int offset = 0; offset < DIGIT_SIZE; ++offset) {
            *(ptr++) = static_cast<int>(dig % 10);
            dig /= 10;
        }
    }
    while (digitsBigEndian.size() > 1 && digitsBigEndian.back() == 0) {
        digitsBigEndian.pop_back();
    }
    std::reverse(digitsBigEndian.begin(), digitsBigEndian.end());
    *this = 0;
    BigInteger cur;
    for (int& dig: digitsBigEndian) {
        cur *= 10;
        *this *= 10;
        cur += BigInteger(dig);
        while (cur >= y) {
            cur -= y;
            ++(*this);
        }
    }
    negative = resultingNegativeness;
    normalize();
    if (f) {
        std::cerr << "division done" << std::endl;
        std::cerr << (double(clock()) - cl) / CLOCKS_PER_SEC << std::endl;
    }
    return *this;
}

BigInteger BigInteger::operator%(const BigInteger& y) const {
    BigInteger res(*this);
    return (res %= y);
}

BigInteger& BigInteger::operator%=(const BigInteger& y) {
    if (negative) {
        if (y.negative) {
            return *this = (-(*this)) % (-y);
        }
        return *this = -((-(*this)) % y);
    }
    if (y.negative) {
        return *this = -((*this) % (-y));
    }
    return *this = *this - ((*this) / y) * y;
}

class Rational {
private:
    BigInteger numerator;
    BigInteger denominator = 1;

    static BigInteger gcd(BigInteger x, BigInteger y) {
        if (x < y)
            x.swap(y);
        int cnt2 = 0;
        while (y) {
            if ((x.digit[0]) % 2 == 0) {
                x /= 2;
                if ((y.digit[0]) % 2 == 0) {
                    y /= 2;
                    ++cnt2;
                }
            } else if ((y.digit[0]) % 2 == 0) {
                y /= 2;
            } else {
               x -= y;
            }
            if (x < y)
                x.swap(y);
        }
        for (int i = 0; i < cnt2; ++i)
            x *= 2;
        return x;
    }

public:
    friend std::istream& operator>>(std::istream& in, Rational& x);

    [[nodiscard]] std::string toString() const;

    void swap(Rational& oth);

    Rational operator-();

    Rational& operator=(const Rational& oth);

    Rational& operator+=(const Rational& oth);
    Rational& operator-=(const Rational& oth);
    Rational& operator*=(const Rational& oth);
    Rational& operator/=(const Rational& oth);

    bool operator<(const Rational& oth) const;
    bool operator>(const Rational& oth) const;
    bool operator!=(const Rational& oth) const;
    bool operator==(const Rational& oth) const;
    bool operator<=(const Rational& oth) const;
    bool operator>=(const Rational& oth) const;

    [[nodiscard]] std::string asDecimal(size_t precision = 0) const;

    Rational() : numerator(0) {}

    Rational(int x) : numerator(x) {}

    Rational(long long x) : numerator(x) {}

    Rational(const Rational& x) = default;

    void normalize() {
        BigInteger x;
        numerator.negative ^= denominator.negative;
        denominator.negative = false;
        if (numerator == 0_bi) {
            denominator = 1_bi;
            numerator.negative = false;
            return;
        }
        x = gcd((numerator.negative ? -numerator : numerator), denominator);
        numerator /= x;
        denominator /= x;

    }

    explicit operator double() const {
        return std::stod(asDecimal(18));
    }
};

Rational Rational::operator-() {
    Rational res(*this);
    res.numerator.negative ^= 1;
    res.normalize();
    return res;
}

void Rational::swap(Rational& oth) {
    numerator.swap(oth.numerator);
    denominator.swap(oth.denominator);
}

Rational& Rational::operator=(const Rational& oth) {
    if (this == &oth)
        return *this;
    Rational tmp = oth;
    this->swap(tmp);
    return *this;
}

bool Rational::operator<(const Rational& oth) const {
    return numerator * oth.denominator < denominator * oth.numerator;
}
bool Rational::operator==(const Rational& oth) const {
    return numerator * oth.denominator == denominator * oth.numerator;
}
bool Rational::operator!=(const Rational& oth) const {
    return numerator * oth.denominator != denominator * oth.numerator;
}
bool Rational::operator>(const Rational& oth) const {
    return numerator * oth.denominator > denominator * oth.numerator;
}
bool Rational::operator<=(const Rational& oth) const {
    return numerator * oth.denominator <= denominator * oth.numerator;
}
bool Rational::operator>=(const Rational& oth) const {
    return numerator * oth.denominator >= denominator * oth.numerator;
}

Rational& Rational::operator*=(const Rational& oth) {
    numerator *= oth.numerator;
    denominator *= oth.denominator;
    normalize();
    return *this;
}
Rational& Rational::operator/=(const Rational& oth) {
    numerator *= oth.denominator;
    denominator *= oth.numerator;
    normalize();
    return *this;
}
Rational& Rational::operator+=(const Rational& oth) {
    numerator *= oth.denominator;
    numerator += denominator * oth.numerator;
    denominator *= oth.denominator;
    normalize();
    return *this;
}
Rational& Rational::operator-=(const Rational& oth) {
    numerator *= oth.denominator;
    numerator -= denominator * oth.numerator;
    denominator *= oth.denominator;
    normalize();
    return *this;
}

Rational operator*(const Rational& x, const Rational& y) {
    Rational res(x);
    return res *= y;
}
Rational operator/(const Rational& x, const Rational& y) {
    Rational res(x);
    return res /= y;
}
Rational operator+(const Rational& x, const Rational& y) {

    Rational res(x);
    return res += y;
}
Rational operator-(const Rational& x, const Rational& y) {
    Rational res(x);
    return res -= y;
}

std::string Rational::toString() const {
    std::string res;
    res += numerator.toString();
    if (denominator != 1_bi) {
        res += '/';
        res += denominator.toString();
    }
    return res;
}
std::string Rational::asDecimal(size_t precision) const {
    if (numerator == 0_bi) {
        std::string res = "0";
        if (precision != 0) {
            res += ".";
        }
        res += std::string(precision, '0');
        return res;
    }
    BigInteger x = numerator;
    for (size_t i = 0; i < precision; ++i) {
        x *= 10_bi;
    }
    x /= denominator;
    std::string res;
    if (x.negative) {
        res += "-";
    }
    x.negative = false;
    std::string s = x.toString();
    if (numerator < denominator) {
        res += "0." + std::string(precision - s.size(), '0');
        res += s;
    } else {
        res += s.substr(0, s.size() - precision) + "." + s.substr(s.size() - precision);
    }
    return res;
}
std::istream& operator>>(std::istream& in, Rational& x) {
    in >> x.numerator;
    x.denominator = 1;
    return in;
}

template<size_t N, size_t potentialDivisor> // Правда ли, что N не делится ни на одно нечетное число >= potentialDivisor. Если делится, то какое среди них минимальное?
struct checkPrime_ {
    static const bool prime = checkPrime_<N,
                                          (potentialDivisor * potentialDivisor > N ? 1 : (N % potentialDivisor == 0 ? 0
                                                                                                                    :
                                                                                          potentialDivisor
                                                                                              + 2))>::prime;
};

template<size_t N>
struct checkPrime_<N, 1> {
    static const bool prime = true;
};

template<size_t N>
struct checkPrime_<N, 0> {
    static const bool prime = false;
};

template<size_t N>
struct isPrime {
    static const bool prime = (N % 2 != 0) && (checkPrime_<N, 3>::prime);
};

template<>
struct isPrime<1> {
    static const bool prime = false;
};

template<>
struct isPrime<2> {
    static const bool prime = true;
};

template<size_t N>
class Residue {
private:
    long long value;
    template<size_t M>
    friend std::ostream& operator<<(std::ostream& out, const Residue<M>& x);
public:
    explicit Residue(int x = 0);
    explicit operator int() const;

    bool operator==(const Residue<N>& y) const;
    bool operator!=(const Residue<N>& y) const;

    Residue<N> operator-() const;
    Residue<N>& operator+=(const Residue<N>& y);
    Residue<N>& operator-=(const Residue<N>& y);
    Residue<N>& operator*=(const Residue<N>& y);
    Residue<N>& operator/=(const Residue<N>& y);

    Residue<N> operator+(const Residue<N>& y) const;
    Residue<N> operator-(const Residue<N>& y) const;
    Residue<N> operator*(const Residue<N>& y) const;
    Residue<N> operator/(const Residue<N>& y) const;

    Residue<N> getInverse() const;
};

template<size_t N>
Residue<N> Residue<N>::operator-() const {
    return Residue<N>(N - value);
}

template<size_t N>
Residue<N>::Residue(int x) {
    value = x % static_cast<int>(N);
    if (value < 0)
        value += N;
}

template<size_t N>
Residue<N>::operator int() const {
    return static_cast<int>(value);
}

template<size_t N>
bool Residue<N>::operator==(const Residue<N>& y) const {
    return (value == y.value);
}

template<size_t N>
bool Residue<N>::operator!=(const Residue<N>& y) const {
    return value != y.value;
}

template<size_t N>
Residue<N>& Residue<N>::operator+=(const Residue<N>& y) {
    value += y.value;
    if (value >= static_cast <long long>(N))
        value -= static_cast <long long>(N);
    return *this;
}

template<size_t N>
Residue<N>& Residue<N>::operator-=(const Residue<N>& y) {
    value -= y.value;
    if (value < 0)
        value += N;
    return *this;
}

template<size_t N>
Residue<N>& Residue<N>::operator*=(const Residue<N>& y) {
    value = (value * y.value) % N;
    return *this;
}

template<size_t N>
Residue<N>& Residue<N>::operator/=(const Residue<N>& y) {
    *this *= y.getInverse();
    return *this;
}

template<size_t N>
Residue<N> Residue<N>::operator+(const Residue<N>& y) const {
    Residue<N> result = *this;
    result += y;
    return result;
}

template<size_t N>
Residue<N> Residue<N>::operator-(const Residue<N>& y) const {
    Residue<N> result = *this;
    result -= y;
    return result;
}

template<size_t N>
Residue<N> Residue<N>::operator*(const Residue<N>& y) const {
    Residue<N> result = *this;
    result *= y;
    return result;
}

template<size_t N>
Residue<N> Residue<N>::operator/(const Residue<N>& y) const {
    Residue<N> result = *this;
    result /= y;
    return result;
}

template<size_t N>
Residue<N> Residue<N>::getInverse() const {
    static_assert(isPrime<N>::prime, "Modulo must be prime to get inverse");
    Residue<N> rev = pow(*this, N - 2);
    return rev;
}

template<size_t N>
Residue<N> pow(Residue<N> x, size_t deg) {
    Residue<N> result = Residue<N>(1);
    while (deg) {
        if (deg & 1)
            result *= x;
        x *= x;
        deg >>= 1;
    }
    return result;
}

template<size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& x) {
    out << x.value;
    return out;
}




//----------------------------------------------------------------------------------------------------------------------

template<size_t N, size_t M, typename Field = Rational>
class Matrix {
private:
    std::vector<std::array<Field, M>> arr;
    std::pair<bool, size_t> gaussForwardStep();

    void swapRows(size_t i, size_t j);
    void addRow(size_t addFrom, size_t addTo, Field k);
    void multiplyRow(size_t i, Field k);
    // Elementary row operations

public:
    static Matrix<N, N, Field> identityMatrix() {
        Matrix<N, N, Field> tmp;
        for (size_t i = 0; i < N; ++i)
            tmp.arr[i][i] = Field(1);
        return tmp;
    }
    

    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field>& y);
    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field>& y);

    Matrix<N, M, Field>& operator*=(const Field& y);
    Matrix<N, M, Field>& operator*=(const Matrix<M, M, Field>& y);

    [[nodiscard]] Field det() const;
    [[nodiscard]] Matrix<M, N, Field> transposed() const;
    [[nodiscard]] Matrix<N, N, Field> inverted() const;
    void invert();
    [[nodiscard]] int rank() const;
    [[nodiscard]] Field trace() const;

    std::array<Field, M> getRow(size_t row) const;
    std::array<Field, N> getColumn(size_t col) const;

    std::array<Field, M>& operator[](size_t ind) {
        return arr[ind];
    }
    

    const std::array<Field, M>& operator[](size_t ind) const {
        return arr[ind];
    }
    

    Matrix() {
        std::array <Field, M> I_LOVE_C_PLUS_PLUS;
        I_LOVE_C_PLUS_PLUS.fill(Field(0));
        arr.resize(N, I_LOVE_C_PLUS_PLUS);
    }

    template<typename T>
    explicit Matrix(const std::vector<std::vector<T>>& a) {
        assert(a.size() == N && a.begin()->size() == M);
        std::array <Field, M> I_LOVE_C_PLUS_PLUS;
        I_LOVE_C_PLUS_PLUS.fill(Field(0));
        arr.resize(N, I_LOVE_C_PLUS_PLUS);
        size_t i = 0;
        for (const auto& row: a) {
            size_t j = 0;
            for (auto el: row)
                arr[i][j++] = Field(el);
            ++i;
        }
    }
    

    template<typename T>
    Matrix(const std::initializer_list<std::vector<T>>& a) {
        assert(a.size() == N && a.begin()->size() == M);
        std::array <Field, M> I_LOVE_C_PLUS_PLUS;
        I_LOVE_C_PLUS_PLUS.fill(Field(0));
        arr.resize(N, I_LOVE_C_PLUS_PLUS);
        size_t i = 0;
        for (const auto& row: a) {
            size_t j = 0;
            for (auto el: row)
                arr[i][j++] = Field(el);
            ++i;
        }
    }

    template<typename T>
    Matrix(const std::initializer_list<std::initializer_list<T>>& a) {
        assert(a.size() == N && a.begin()->size() == M);
        std::array <Field, M> I_LOVE_COPY_PASTA;
        I_LOVE_COPY_PASTA.fill(Field(0));
        arr.resize(N, I_LOVE_COPY_PASTA);
        size_t i = 0;
        for (const auto& row: a) {
            size_t j = 0;
            for (auto el: row)
                arr[i][j++] = Field(el);
            ++i;
        }
    }

    void print() {
        for (size_t i = 0; i < N; ++i) {
            for (const auto& el: getRow(i))
                cout << el.toString() << ' ';
            cout << endl;
        }
    }
    
};

template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator+=(const Matrix<N, M, Field>& y) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            arr[i][j] += y.arr[i][j];
        }
    }
    return *this;
}


template<size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator-=(const Matrix<N, M, Field>& y) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            arr[i][j] -= y.arr[i][j];
        }
    }
    return *this;
}


template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    Matrix<N, M, Field> ans = x;
    ans += y;
    return ans;
}


template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    Matrix<N, M, Field> ans = x;
    ans -= y;
    return ans;
}


template<size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(const Field& y) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            arr[i][j] *= y;
        }
    }
    return *this;
}


template<size_t N, size_t M, typename Field>
Matrix<M, N, Field> Matrix<N, M, Field>::transposed() const {
    Matrix<M, N, Field> ans;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            ans[i][j] = arr[j][i];
        }
    }
    return ans;
}


template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& x, const Field& k) {
    Matrix<N, M, Field> ans = x;
    ans *= k;
    return ans;
}


template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Field& k, const Matrix<N, M, Field>& x) {
    return x * k;
}


template<size_t N, size_t M, size_t L, typename Field>
Matrix<N, L, Field> operator*(const Matrix<N, M, Field>& x, const Matrix<M, L, Field>& y) {
    Matrix<N, L, Field> ans;
    for (size_t z = 0; z < M; ++z) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < L; ++j) {
                ans[i][j] += x[i][z] * y[z][j];
            }
        }
    }
    return ans;
}


template<size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(const Matrix<M, M, Field>& y) {
    static_assert(N == M, "Matrices must be square to perform *=");
    return *this = *this * y;
}


template<size_t N, size_t M, typename Field>
std::array<Field, M> Matrix<N, M, Field>::getRow(size_t row) const {
    return arr[row];
}


template<size_t N, size_t M, typename Field>
std::array<Field, N> Matrix<N, M, Field>::getColumn(size_t col) const {
    std::array <Field, N> ans;
    for (size_t i = 0; i < N; ++i)
        ans[i] = arr[i][col];
    return ans;
}


template<size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::trace() const {
    static_assert(N == M, "Matrix must be square to calculate its trace");
    Field ans(0);
    for (size_t i = 0; i < N; ++i) {
        ans += arr[i][i];
    }
    return ans;
}


template<size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::det() const {
    static_assert(N == M, "Matrix must be square to calculate its determinant");
    Matrix<N, M, Field> copy = *this;
    bool inv = copy.gaussForwardStep().first;
    Field result = inv ? Field(-1) : Field(1);
    for (size_t i = 0; i < N; ++i) {
        result *= copy.arr[i][i];
    }
    return result;
}


template<size_t N, size_t M, typename Field>
Matrix<N, N, Field> Matrix<N, M, Field>::inverted() const {
    Matrix<N, M, Field> copy = *this;
    copy.invert();
    return copy;
}


template<size_t N, size_t M, typename Field>
void Matrix<N, M, Field>::invert() {
    static_assert(N == M, "Matrix must be square to invert it");
    // Приставляем нас к единичной матрице ans
    auto ans = Matrix<N, N, Field>::identityMatrix();

    for (size_t row = 0, col = 0; row < N && col < M; ++col) {
        bool isColumnZero = false;
        if (arr[row][col] == Field(0)) {
            isColumnZero = true;
            for (size_t rowToSwap = row + 1; rowToSwap < N; ++rowToSwap) {
                if (arr[rowToSwap][col] != Field(0)) {
                    this->swapRows(row, rowToSwap);
                    ans.swapRows(row, rowToSwap);
                    isColumnZero = false;
                    break;
                }
            }
        }
        assert(!isColumnZero);
        // Инвертируемая матрица обязана быть невырожденной
        for (size_t followingRow = row + 1; followingRow < N; ++followingRow) {
            if (arr[followingRow][col] != Field(0)) {
                Field k = (-arr[followingRow][col] / arr[row][col]);
                this->addRow(row, followingRow, k);
                ans.addRow(row, followingRow, k);
            }
        }
        ++row;
    }
    //gauss backward step
    for (int col = N - 1; col >= 0; --col) {
        int& row = col;
        Field k = Field(1) / arr[row][col];
        this->multiplyRow(row, k);
        ans.multiplyRow(row, k);
        for (int precedingRow = row - 1; precedingRow >= 0; --precedingRow) {
            if (arr[precedingRow][col] == Field(0))
                continue;
            k = (-arr[precedingRow][col]);
            this->addRow(row, precedingRow, k);
            ans.addRow(row, precedingRow, k);
        }
    }
    *this = ans;
}


template<size_t N, size_t M, typename Field>
int Matrix<N, M, Field>::rank() const {
    Matrix<N, M, Field> copy = *this;
    return copy.gaussForwardStep().second;
}

template<size_t N, size_t M, typename Field>
void Matrix<N, M, Field>::swapRows(size_t i, size_t j) {
    arr[i].swap(arr[j]);
}


template<size_t N, size_t M, typename Field>
void Matrix<N, M, Field>::addRow(size_t addFrom, size_t addTo, Field k) {
    for (size_t col = 0; col < M; ++col) {
        arr[addTo][col] += arr[addFrom][col] * k;
    }
}


template<size_t N, size_t M, typename Field>
void Matrix<N, M, Field>::multiplyRow(size_t i, Field k) {
    for (size_t j = 0; j < M; ++j) {
        arr[i][j] *= k;
    }
}


// returns parity of amount of swaps done and rank of matrix
template<size_t N, size_t M, typename Field>
std::pair<bool, size_t> Matrix<N, M, Field>::gaussForwardStep() {
    size_t rank_ = 0;
    bool parityOfSwaps = false;
    size_t row = 0;
    size_t col = 0;
    for (; row < N && col < M; ++col) {
        bool isColumnZero = false;
        if (arr[row][col] == Field(0)) {
            isColumnZero = true;
            for (size_t rowToSwap = row + 1; rowToSwap < N; ++rowToSwap) {
                if (arr[rowToSwap][col] != Field(0)) {
                    swapRows(row, rowToSwap);
                    parityOfSwaps ^= 1;
                    isColumnZero = false;
                    break;
                }
            }
        }
        if (isColumnZero)
            continue;
        for (size_t followingRow = row + 1; followingRow < N; ++followingRow) {
            if (arr[followingRow][col] != Field(0))
                addRow(row, followingRow, -arr[followingRow][col] / arr[row][col]);
        }
        ++rank_; // Rank equals to amount of row echelons
        ++row;
    }
    return {parityOfSwaps, rank_};
}


template<size_t N, size_t M, typename Field>
bool operator==(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            if (x[i][j] != y[i][j]) {
                return false;
            }
        }
    }
    return true;
}


template<size_t N, size_t M, typename Field>
bool operator!=(const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
    return !(x == y);
}



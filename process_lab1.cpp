// Лабораторная работа №1. Вариант 1б.
// Многократные прямые измерения физических величин
// и обработка результатов наблюдений.
//
// Кафедра физической механики, СПбГУ, 2026
// Компиляция: g++ -std=c++17 -o process_lab1 process_lab1.cpp

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <numeric>

// ─── Функции статистической обработки ────────────────────────────────────────

// Среднее арифметическое (формула 5 методички)
double mean(const std::vector<double>& x)
{
    double s = 0.0;
    for (double v : x) s += v;
    return s / static_cast<double>(x.size());
}

// Отклонения d_i = x_i - x_mean
std::vector<double> deviations(const std::vector<double>& x, double x_mean)
{
    std::vector<double> d(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        d[i] = x[i] - x_mean;
    return d;
}

// Сумма квадратов отклонений
double sum_sq(const std::vector<double>& d)
{
    double s = 0.0;
    for (double v : d) s += v * v;
    return s;
}

// Стандартная погрешность отдельного наблюдения (формула 14)
double std_dev(const std::vector<double>& d)
{
    int n = static_cast<int>(d.size());
    return std::sqrt(sum_sq(d) / (n - 1));
}

// Погрешность среднеарифметического (формула 15)
double std_mean(double s, int n)
{
    return s / std::sqrt(static_cast<double>(n));
}

// Коэффициент Стьюдента: таблица для alpha=0.95
// (значения из приложения 3 методического пособия)
double student_t(int n)
{
    // n — число измерений; ищем t по числу степеней свободы f = n-1
    // Таблица: f -> t(0.95)
    const int   fs[] = { 1,  2,    3,    4,    5,    6,    7,    8,
                         9, 10,   15,   20,   25,   30,   40,   50, 1000};
    const double ts[] = {12.71, 4.30, 3.18, 2.78, 2.57, 2.45, 2.36, 2.31,
                          2.26, 2.23, 2.13, 2.09, 2.06, 2.04, 2.02, 2.01, 1.96};
    int f = n - 1;
    // Поиск ближайшего значения
    int idx = 0;
    for (int i = 0; i < 17; ++i) {
        if (f >= fs[i]) idx = i;
        else break;
    }
    return ts[idx];
}

// Доверительный интервал (формула 17)
double conf_interval(double t, double s_mean_val)
{
    return t * s_mean_val;
}

// Погрешность прибора вольтметра (формула из описания прибора)
// S = ±(a + b * Uk/Ux)%,  omega = S/100 * Ux
double instrument_error(double Uk, double Ux, double a = 0.05, double b = 0.05)
{
    double S_pct = a + b * (Uk / Ux);
    return S_pct / 100.0 * Ux;
}

// Суммарная погрешность (ситуация 3, формула 18)
double total_error(double omega, double s_mean_val)
{
    double w3 = omega / 3.0;
    return std::sqrt(w3 * w3 + s_mean_val * s_mean_val);
}

// ─── Гистограмма ─────────────────────────────────────────────────────────────

struct Bin {
    double low, high, mid;
    int    count;
    double freq;  // delta_n / n
};

std::vector<Bin> histogram(const std::vector<double>& x,
                           double lo, double hi, int num_bins)
{
    int n = static_cast<int>(x.size());
    double width = (hi - lo) / num_bins;
    std::vector<Bin> bins(num_bins);
    for (int k = 0; k < num_bins; ++k) {
        bins[k].low   = lo + k * width;
        bins[k].high  = lo + (k + 1) * width;
        bins[k].mid   = (bins[k].low + bins[k].high) / 2.0;
        bins[k].count = 0;
    }
    for (double v : x) {
        int k = static_cast<int>((v - lo) / width);
        if (k >= 0 && k < num_bins) bins[k].count++;
        else if (v == hi)            bins[num_bins - 1].count++;
    }
    for (auto& b : bins)
        b.freq = static_cast<double>(b.count) / n;
    return bins;
}

// ─── Вывод данных для gnuplot ─────────────────────────────────────────────────

void write_timeseries(const std::string& fname,
                      const std::vector<double>& x, double x_mean)
{
    std::ofstream f(fname);
    f << "# i  U_i  U_mean\n";
    for (size_t i = 0; i < x.size(); ++i)
        f << (i + 1) << "\t" << std::fixed << std::setprecision(4)
          << x[i] << "\t" << x_mean << "\n";
}

void write_histogram(const std::string& fname,
                     const std::vector<Bin>& bins)
{
    std::ofstream f(fname);
    f << "# mid  delta_n/n\n";
    for (const auto& b : bins)
        f << std::fixed << std::setprecision(4) << b.mid << "\t"
          << std::setprecision(4) << b.freq << "\n";
}

// ─── Главная программа ────────────────────────────────────────────────────────

int main()
{
    std::cout << std::fixed;

    // ── Данные ──────────────────────────────────────────────────────────────
    const std::vector<double> data = {
        0.3506, 0.3508, 0.3516, 0.3509, 0.3507,
        0.3506, 0.3505, 0.3489, 0.3498, 0.3500,
        0.3503, 0.3504, 0.3501, 0.3503, 0.3507,
        0.3508, 0.3507, 0.3503, 0.3508, 0.3505,
        0.3506, 0.3508, 0.3509, 0.3508, 0.3509,
        0.3510, 0.3511, 0.3509, 0.3508, 0.3514,
        0.3516, 0.3511, 0.3500, 0.3499, 0.3501,
        0.3503, 0.3504, 0.3503, 0.3497, 0.3501,
        0.3511, 0.3508, 0.3509, 0.3510, 0.3511,
        0.3508, 0.3514, 0.3511, 0.3512, 0.3506
    };
    const int n = static_cast<int>(data.size());

    // ── Параметры прибора ────────────────────────────────────────────────────
    const double Uk_fine   =  1.0;    // точная шкала, В

    // ── Статистическая обработка ─────────────────────────────────────────────
    double x_mean = mean(data);
    auto   d      = deviations(data, x_mean);
    double ss     = sum_sq(d);
    double s      = std_dev(d);
    double sm     = std_mean(s, n);
    double t      = student_t(n);
    double eps    = conf_interval(t, sm);
    double omega  = instrument_error(Uk_fine, x_mean);
    double omega3 = omega / 3.0;
    double s_sum  = total_error(omega, sm);

    // ── Вывод таблицы d_i ───────────────────────────────────────────────────
    std::cout << "=== ТАБЛИЦА НАБЛЮДЕНИЙ (точная шкала 1 В) ===\n";
    std::cout << std::setw(5)  << "i"
              << std::setw(10) << "U_i, В"
              << std::setw(14) << "d_i, В"
              << std::setw(16) << "d_i^2, В^2"  << "\n";
    std::cout << std::string(45, '-') << "\n";
    for (int i = 0; i < n; ++i) {
        std::cout << std::setw(5)  << (i + 1)
                  << std::setw(10) << std::setprecision(4) << data[i]
                  << std::setprecision(6)
                  << std::setw(14) << d[i]
                  << std::setprecision(8)
                  << std::setw(16) << d[i] * d[i] << "\n";
    }
    std::cout << std::string(45, '-') << "\n";
    std::cout << "Сумма d_i^2 = " << std::setprecision(8) << ss << " В^2\n\n";

    // ── Вывод результатов расчёта ────────────────────────────────────────────
    std::cout << "=== РЕЗУЛЬТАТЫ РАСЧЁТА ===\n";
    std::cout << "n                  = " << n << "\n";
    std::cout << "U_mean             = " << std::setprecision(6) << x_mean << " В\n";
    std::cout << "s (sigma_single)   = " << std::setprecision(6) << s      << " В\n";
    std::cout << "s_mean             = " << std::setprecision(6) << sm     << " В\n";
    std::cout << "t(alpha=0.95,n=50) = " << std::setprecision(2) << t      << "\n";
    std::cout << "eps (delta_rand)   = " << std::setprecision(6) << eps    << " В\n";
    std::cout << "omega (delta_instr)= " << std::setprecision(6) << omega  << " В\n";
    std::cout << "omega/3            = " << std::setprecision(6) << omega3 << " В\n";
    std::cout << "s / (omega/3)      = " << std::setprecision(3) << s / omega3 << "\n";
    std::cout << "Ситуация           = 3  (omega/3 < s < omega*4 — оба вклада)\n";
    std::cout << "s_sum              = " << std::setprecision(6) << s_sum  << " В\n\n";
    std::cout << "ОКОНЧАТЕЛЬНЫЙ РЕЗУЛЬТАТ (alpha = 0.95):\n";
    std::cout << "  U = (" << std::setprecision(4) << x_mean
              << " +/- " << std::setprecision(4) << eps << ") В\n\n";

    // ── Гистограмма ─────────────────────────────────────────────────────────
    auto bins = histogram(data, 0.3488, 0.3520, 8);
    std::cout << "=== ГИСТОГРАММА ===\n";
    std::cout << std::setw(20) << "Интервал, В"
              << std::setw(8)  << "dn"
              << std::setw(12) << "dn/n" << "\n";
    std::cout << std::string(40, '-') << "\n";
    for (const auto& b : bins) {
        std::cout << "[" << std::setprecision(4) << b.low
                  << "; " << b.high << ")"
                  << std::setw(6) << b.count
                  << std::setw(12) << std::setprecision(4) << b.freq << "\n";
    }

    // ── Данные для gnuplot ───────────────────────────────────────────────────
    write_timeseries("timeseries.dat", data, x_mean);
    write_histogram("histogram.dat", bins);
    std::cout << "\nФайлы timeseries.dat и histogram.dat записаны.\n";

    return 0;
}

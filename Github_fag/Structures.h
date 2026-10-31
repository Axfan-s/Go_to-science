#pragma once
#include <vector>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <stdexcept>
#include <conio.h>
#define NOMINMAX
#include <windows.h>
#include <queue>
//#include "Mutations.h"
// Исправленная генерация случайных чисел
std::vector<long double> Zero_Copy(std::vector <long double> Cur, int k) {
    std::vector<long double> New;
    for (int i = 0; i < Cur.size(); ++i) {
        if (i != k) {
            New.push_back(Cur[i]);
        }
        else {
            New.push_back(0);
        }
    }
    return New;
}
void saveTwoVectorsToCSV(const std::vector<double>& vec1,
    const std::vector<double>& vec2,
    const std::string& filename,
    const std::string& col1_name = "Column1",
    const std::string& col2_name = "Column2") {

    std::ofstream file(filename);

    file << "Vector1" << ";" << "Vector2\n";

    size_t min_size = std::min(vec1.size(), vec2.size());
    for (size_t i = 0; i < min_size; ++i) {
        file << vec1[i] << ";" << vec2[i] << "\n";
    }

    file.close();
}

void saveSetVectorsToCSV(std::vector<std::vector<long double>> data,
    const std::string& filename, const std::string& col_name = "Column") {

    std::ofstream file(filename);
    for (int i = 0; i < data.size(); ++i) {
        if (i != data.size() - 1) {
            file << col_name << i + 1 << ";";
        }
        else {
            file << col_name << i + 1 << "\n";
        }
    }
    size_t sizet = data[0].size();
    for (int j = 0; j < sizet; ++j) {

        for (size_t i = 0; i < data.size(); ++i) {
            if (i != data.size() - 1) {
                file << "\u200B" << data[i][j] << ";";
            }
            else {
                file << "\u200B" << data[i][j] << "\n";
            }
        }
    }
    //file.close();
}
double my_random_d(double a = 1, double b = 0) {
    static std::random_device rd;
    static std::mt19937_64 generator(rd());

    if (b > a) std::swap(a, b);
    std::uniform_real_distribution<double> dist(b, a);
    return dist(generator);
}
int my_random_i(int a = 1000, int b = 0) {
    static std::random_device rd;
    static std::mt19937_64 generator(rd());

    if (b > a) std::swap(a, b);
    std::uniform_int_distribution<int> dist(b, a);
    return dist(generator);
}
unsigned long seed = 1;

/*void init_rand() {
    srand(time(0));
    seed = rand();
}

int my_rand() {
    seed = seed * 1103515245 + 12345;
    return (seed / 655) % 1000000;
}*/

class Point {
private:
    size_t dimension;
    std::vector< long double> point;

public:

    // Конструкторы
    Point() : dimension(0), point(0, 0.0) {}

    Point(size_t dim) : dimension(dim), point(dim, 0.0) {}

    Point(const std::vector<long double>& coords) : dimension(coords.size()), point(coords) {}

    operator std::vector<long double>() const {
        return point;
    }
    Point operator+(const Point& A) const {
        if (this->dimension != A.dimension) {
            return Point();
        }

        Point result(this->dimension);
        for (size_t i = 0; i < dimension; ++i) {
            result[i] = this->point[i] + A.point[i];
        }
        return result;
    }
    Point operator+(int A) const {
        Point result(this->dimension);
        for (size_t i = 0; i < dimension; ++i) {
            result[i] = this->point[i] + A;
        }
        return result;
    }
    Point operator-(const Point& A) const {
        if (this->dimension != A.dimension) {
            return Point();
        }

        Point result(this->dimension);
        for (size_t i = 0; i < dimension; ++i) {
            result[i] = this->point[i] - A.point[i];
        }
        return result;
    }

    Point operator*(double a) const {
        Point result(this->dimension);
        for (size_t i = 0; i < dimension; ++i) {
            result[i] = a * this->point[i];
        }
        return result;
    }

    Point operator~() {
        Point result(dimension);
        for (size_t i = 0; i < dimension; ++i) {
            result[i] = !point[i];
        }
        return result;
    }
    // Оператор индексирования
    long double& operator[](size_t index) {
        return point[index];
    }

    const long double& operator[](size_t index) const {
        return point[index];
    }

    // Метод для вычисления нормы
    long double norm() const {
        long double square = 0;
        for (size_t i = 0; i < dimension; ++i) {
            square += point[i] * point[i];
        }
        return sqrt(square);
    }
    void Normalization() {
        long double len = this->norm();  // или просто norm()
        if (len > 0) {  // защита от деления на ноль
            for (size_t i = 0; i < dimension; ++i) {
                point[i] = point[i] / len;
            }
        }
    }

    // Метод для генерации точки
    void MakePoint(double a, double b) {
        bool l = true;
        for (size_t i = 0; i < dimension; ++i) {
            double R = my_random_d(a, b);
            point[i] = R;
        }
    }

    // Геттеры
    size_t getDimension() const { return dimension; }
    const std::vector<long double>& getPoint() const { return point; }

    /* void generate_mutation(const DistributionParams& params) {
         switch (params.type) {
         case DistributionParams::NORMAL:
             for (size_t i = 0; i < this->getDimension(); ++i) {
                 double q = my_rand() % 1001 * 1.0 / 1000;
                 double eta = my_rand() % 1001 * 1.0 / 1000;
                 double Expr = eta * (boost::math::erf((*params.param2 - (*this)[i]) / sqrt(2 * *params.param3))) + (1 - eta) * boost::math::erf((*params.param1 - (*this)[i]) / sqrt(2 * *params.param3));
                 double R;
                 if (Expr >= 1) {
                     R = *params.param2;
                 }
                 else if (Expr <= -1) {
                     R = *params.param1;
                 }
                 else {
                     R = sqrt(2 * *params.param3) * (boost::math::erf_inv(Expr)) + (*this)[i];
                 }
                 (*this)[i] = R;
             }
             return;
             break;

         case DistributionParams::CAUCHY:
             if (params.param1 && params.param2) {
                 std::cauchy_distribution<double> dist(*params.param1, *params.param2);
                 return ;
             }
             break;

         case DistributionParams::UNIFORM:
             if (params.param1 && params.param2) {
                 std::uniform_real_distribution<double> dist(*params.param1, *params.param2);
                 return ;
             }
             break;

         case DistributionParams::EXPONENTIAL:
             if (params.param1) {
                 std::exponential_distribution<double> dist(*params.param1);
                 return ;
             }
             break;

         case DistributionParams::LAPLACE:
             if (params.param1 && params.param2) {
                 // Реализация распределения Лапласа
                 std::uniform_real_distribution<double> uniform(0.0, 1.0);
                 double u = 0.5;
                 return ;
             }
             break;
         }
         return ;
     }*/
     /*void mutate(const DistributionParams& params) {
         for (size_t i = 0; i < dimension; ++i) {
             point[i] = generate_mutation(params);
         }
     }*/
    long double distanceTo(Point A) {
        Point C = A - *this;
        return C.norm();
    }
};

long double distance(Point A, Point B) {
    Point C = A - B;
    return C.norm();
}
long double ScalarMultiply(Point a, Point b) {
    long double result = 0;
    for (int i = 0; i < a.getDimension(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}


class PointSet : public Point {
private:
    size_t index;                                    // собственный индекс
    static constexpr size_t DEFAULT_K = 20;          // храним 20 ближайших

    // Храним пары (расстояние, индекс соседа)
    using Neighbor = std::pair<long double, size_t>;
    std::priority_queue<Neighbor> nearest;           // max-heap

public:
    // Конструкторы
    PointSet() : Point(), index(0) {}

    PointSet(size_t dim, size_t idx) : Point(dim), index(idx) {}

    PointSet(const std::vector<long double>& coords, size_t idx)
        : Point(coords), index(idx) {
    }

    // Копирование
    PointSet(const PointSet& other)
        : Point(other), index(other.index), nearest(other.nearest) {
    }

    // Геттеры
    size_t getIndex() const { return index; }

    // Добавление соседа с поддержанием k ближайших
    void addNeighborIfCloser(size_t neighbor_idx, long double distance, size_t k = DEFAULT_K) {
        if (nearest.size() < k) {
            nearest.emplace(distance, neighbor_idx);
        }
        else if (distance < nearest.top().first) {
            nearest.pop();
            nearest.emplace(distance, neighbor_idx);
        }
    }

    // Получение k ближайших соседей
    std::vector<size_t> getKNearest(size_t k = DEFAULT_K) const {
        std::vector<Neighbor> temp;
        temp.reserve(nearest.size());

        // Копируем и извлекаем
        auto copy = nearest;
        while (!copy.empty()) {
            temp.push_back(copy.top());
            copy.pop();
        }

        // Сортируем по возрастанию расстояния
        std::sort(temp.begin(), temp.end());

        k = std::min(k, temp.size());
        std::vector<size_t> result;
        result.reserve(k);
        for (size_t i = 0; i < k; ++i) {
            result.push_back(temp[i].second);
        }
        return result;
    }

    // Получить всех соседей (для отладки)
    std::vector<size_t> getAllNeighbors() const {
        return getKNearest(nearest.size());
    }

    // Очистка соседей
    void clearNeighbors() {
        nearest = std::priority_queue<Neighbor>();
    }

    // Проверка, есть ли соседи
    bool hasNeighbors() const { return !nearest.empty(); }

    // Количество хранимых соседей
    size_t neighborCount() const { return nearest.size(); }
};

class Population {
private:
    Point* points;
    size_t size;
    size_t dimension;

public:
    Population() : points(nullptr), size(0), dimension(0) {}

    Population(size_t size, size_t dimension, Point* points)
        : size(size), dimension(dimension) {
        this->points = new Point[size];
        for (size_t i = 0; i < size; ++i) {
            this->points[i] = points[i];
        }
    }

    Population(size_t size, size_t dimension)
        : size(size), dimension(dimension) {
        this->points = new Point[size]; // ВЫДЕЛЯЕМ ПАМЯТЬ!
        for (size_t i = 0; i < size; ++i) {
            this->points[i] = Point(dimension); // Создаем точки нужной размерности
        }
    }

    // Конструктор копирования
    Population(const Population& other)
        : size(other.size), dimension(other.dimension) {
        points = new Point[size];
        for (size_t i = 0; i < size; ++i) {
            points[i] = other.points[i];
        }
    }

    // Оператор присваивания
    Population& operator=(const Population& other) {
        if (this != &other) {
            delete[] points; // Освобождаем старую память

            size = other.size;
            dimension = other.dimension;
            points = new Point[size];
            for (size_t i = 0; i < size; ++i) {
                points[i] = other.points[i];
            }
        }
        return *this;
    }

    // Деструктор
    ~Population() {
        delete[] points; // ОСВОБОЖДАЕМ ПАМЯТЬ!
    }

    // Оператор индексирования
    Point& operator[](size_t index) {
        if (index >= size) {
            throw std::out_of_range("Index out of range");
        }
        return points[index];
    }

    const Point& operator[](size_t index) const {
        if (index >= size) {
            throw std::out_of_range("Index out of range");
        }
        return points[index];
    }

    void GeneratePopulation(double a, double b) {
        for (size_t i = 0; i < size; ++i) {
            points[i].MakePoint(a, b);
        }
    }
    std::vector<std::vector<long double>> PopulationToMatrix() const {
        size_t rows = this->getSize();
        size_t cols = this->getDimension();

        std::vector<std::vector<long double>> Matrix;
        Matrix.reserve(rows);

        for (size_t i = 0; i < rows; ++i) {
            std::vector<long double> row;
            row.reserve(cols);
            for (size_t j = 0; j < cols; ++j) {
                row.push_back((*this)[i][j]);
            }
            Matrix.push_back(std::move(row));
        }

        return Matrix;
    }
    // Геттеры
    size_t getSize() const { return size; }
    size_t getDimension() const { return dimension; }

    int BestAgent(long double (*J)(Point)) {
        double min = 10000000000.0;
        int min_index = INT_MAX;
        double F;
        for (int i = 0; i < size; ++i) {
            F = J(points[i]);
            if (F < min) {
                min = F;
                min_index = i;
            }
        }
        return min_index;
    }

};

void RecordResult(std::vector<Point> Total) {
    std::ofstream file("example.csv");

    if (!file.is_open()) {
        std::cerr << "Failed to open file!" << std::endl;
    }
    for (int i = 0; i < Total.size(); ++i) {
        for (size_t j = 0; j < Total[i].getDimension(); ++j) {
            file << Total[i][j];
            if (i != Total[i].getDimension() - 1) file << ",";
        }
        file << "\n";
    }

    file.close();
    std::cout << "CSV file created successfully." << std::endl;
}

class PopulationSet {
private:
    std::vector<PointSet> points;
    size_t dimension;
    size_t default_k;                    // число соседей по умолчанию

public:
    // Конструктор
    PopulationSet(size_t N, size_t dim, size_t k = 10)
        : dimension(dim), default_k(k) {
        points.reserve(N);
        for (size_t i = 0; i < N; ++i) {
            points.emplace_back(dim, i);  // используем существующий конструктор PointSet(dim, idx)
        }
    }

    // Доступ к точкам
    PointSet& operator[](size_t i) { return points[i]; }
    const PointSet& operator[](size_t i) const { return points[i]; }

    // Размер популяции
    size_t size() const { return points.size(); }

    // Размерность
    size_t getDimension() const { return dimension; }

    // Генерация случайной популяции
    void generatePopulation(double a, double b) {
        for (auto& p : points) {
            p.MakePoint(a, b);           // унаследован от Point
        }
        rebuildNeighbors();
    }

    // Полное перестроение списков соседей (всех точек)
    void rebuildNeighbors() {
        size_t N = points.size();

        // Очищаем всех
        for (auto& p : points) {
            p.clearNeighbors();          // есть в PointSet
        }

        // Вычисляем расстояния
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                if (i == j) continue;
                // distanceTo унаследован от Point
                long double dist = points[i].distanceTo(points[j]);
                points[i].addNeighborIfCloser(j, dist, default_k);  // есть в PointSet
            }
        }
    }

    // Обновление соседей для одной точки (после изменения её координат)
    void updateNeighborsForPoint(size_t idx) {
        points[idx].clearNeighbors();

        for (size_t j = 0; j < points.size(); ++j) {
            if (idx == j) continue;
            long double dist = points[idx].distanceTo(points[j]);
            points[idx].addNeighborIfCloser(j, dist, default_k);
        }
    }

    // Обновление всех точек после изменения точки idx
    void updateAllAfterPointChange(size_t idx) {
        // Обновляем изменённую точку
        updateNeighborsForPoint(idx);

        // Обновляем все остальные точки (у них изменилось расстояние до idx)
        for (size_t i = 0; i < points.size(); ++i) {
            if (i == idx) continue;
            // Проще перестроить для каждой точки заново
            points[i].clearNeighbors();
            for (size_t j = 0; j < points.size(); ++j) {
                if (i == j) continue;
                long double dist = points[i].distanceTo(points[j]);
                points[i].addNeighborIfCloser(j, dist, default_k);
            }
        }
    }

    // Вычисление градиента в точке idx через k ближайших соседей
    Point computeGradient(size_t idx, size_t k, const std::vector<long double>& fitness) {
        Point grad(dimension);

        if (k == 0) return grad;

        auto neighbors = points[idx].getKNearest(k);  // есть в PointSet
        if (neighbors.empty()) return grad;

        for (size_t neighbor_idx : neighbors) {
            long double df = fitness[neighbor_idx] - fitness[idx];
            long double dist = points[idx].distanceTo(points[neighbor_idx]);

            if (std::abs(dist) < 1e-12) continue;

            for (size_t d = 0; d < dimension; ++d) {
                long double diff = points[neighbor_idx][d] - points[idx][d];
                grad[d] += (df / dist) * (diff / dist);
            }
        }

        return grad;
    }

    // Адаптивный градиент (k = min(запрошенное, доступное))
    Point computeAdaptiveGradient(size_t idx, size_t desired_k, const std::vector<long double>& fitness) {
        size_t available = points[idx].neighborCount();  // есть в PointSet
        size_t k = std::min(desired_k, available);
        return computeGradient(idx, k, fitness);
    }

    // Получить индексы k ближайших соседей для точки idx
    std::vector<size_t> getKNearest(size_t idx, size_t k) const {
        return points[idx].getKNearest(k);
    }

    // Получить расстояния до k ближайших соседей (нужно добавить метод в PointSet)
    // Если в PointSet нет метода getKNearestDistances, этот метод нужно убрать
    std::vector<long double> getKNearestDistances(size_t idx, size_t k) const {
        ;
        return {};
    }

    // Вывод всей популяции (для отладки)
    void print(std::ostream& os = std::cout) const {
        os << "Population: " << size() << " points, dimension " << dimension
            << ", neighbors per point: " << default_k << std::endl;
        for (size_t i = 0; i < points.size(); ++i) {
            os << "  Point " << points[i].getIndex() << ": neighbors: "
                << points[i].neighborCount() << "/" << default_k << std::endl;
        }
    }
};

class Domain {
private:
    int dim;
    double radius;
    Point local_min;
    long double (*objective_Func)(Point);
public:
    bool InDomain(Point z) {
        long double dist = distance(z, local_min);
        if (dist > radius) {
            return false;
        }
        else {
            return true;
        }
    }
};

class Menu {
private:
    std::vector<std::string> items;
    int selected = 0;
    std::string title;

    void gotoxy(int x, int y) {
        COORD coord;
        coord.X = x;
        coord.Y = y;
        SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord);
    }

    void setColor(int color) {
        SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
    }

public:
    Menu(const std::string& title) : title(title) {}

    void addItem(const std::string& item) {
        items.push_back(item);
    }

    int run() {
        int key;
        bool running = true;

        while (running) {
            system("cls"); // очистка экрана

            // Заголовок
            std::cout << title << std::endl;

            // Вывод пунктов меню
            for (size_t i = 0; i < items.size(); i++) {
                if (i == selected) {
                    setColor(112); // выделенный пункт (ярко-белый на синем)
                    std::cout << " > " << items[i] << " <" << std::endl;
                    setColor(7);   // сброс цвета
                }
                else {
                    std::cout << "   " << items[i] << std::endl;
                }
            }

            key = _getch();

            if (key == 224) {
                key = _getch();
                switch (key) {
                case 72:
                    if (selected > 0) selected--;
                    break;
                case 80:
                    if (selected < items.size() - 1) selected++;
                    break;
                }
            }
            else if (key == 13) { // Enter
                running = false;
            }
            else if (key == 27) { // Escape
                return -1; // отмена
            }
        }
        return selected;
    }
};



std::vector<long double> solveLinearSystem(const std::vector<Point>& A,
    const std::vector<long double>& b) {
    int n = A.size();

    // Проверка размерности
    if (n == 0 || A[0].getDimension() != n || b.size() != n) {
        throw std::runtime_error("Invalid matrix dimensions");
    }

    // Копируем матрицу и вектор (будем их изменять)
    std::vector<Point> M = A;
    std::vector<long double> v = b;

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; ++i) {
        // Поиск главного элемента (для устойчивости)
        int maxRow = i;
        long double maxVal = std::abs(M[i][i]);
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(M[k][i]) > maxVal) {
                maxVal = std::abs(M[k][i]);
                maxRow = k;
            }
        }

        // Проверка на вырожденность
        /*if (maxVal < 1e-12) {
            throw std::runtime_error("Singular matrix");
        }*/

        // Перестановка строк
        if (maxRow != i) {
            std::swap(M[i], M[maxRow]);
            std::swap(v[i], v[maxRow]);
        }

        // Нормализация текущей строки
        long double pivot = M[i][i];
        for (int j = i; j < n; ++j) {
            M[i][j] /= pivot;
        }
        v[i] /= pivot;

        // Исключение переменной из остальных строк
        for (int k = 0; k < n; ++k) {
            if (k != i && std::abs(M[k][i]) > 1e-12) {
                long double factor = M[k][i];
                for (int j = i; j < n; ++j) {
                    M[k][j] -= factor * M[i][j];
                }
                v[k] -= factor * v[i];
            }
        }
    }

    return v;  // v теперь содержит решение
}


// Структура для хранения расстояния и индекса
typedef struct DistIndex {
    size_t index;
    long double distance;

    DistIndex() : index(0), distance(0.0) {}
    DistIndex(size_t idx, long double dist) : index(idx), distance(dist) {}

    size_t getIndex() const { return index; }
    long double getDist() const { return distance; }
};

// Функция сравнения для сортировки
bool CompareByDist(const DistIndex& A, const DistIndex& B) {
    return A.getDist() < B.getDist();
}
std::vector<DistIndex> getKNearest(const std::vector<DistIndex>& D, size_t k) {
    //volatile bool stop_here = true;
    k = std::min(k, D.size());
    std::vector<DistIndex> result;
    result.reserve(k);
    for (size_t i = 0; i < k; ++i) {
        result.push_back(D[i]);
    }
    return result;
}
std::vector<Point> MatrixToSolve(const Point& point,
    const std::vector<DistIndex>& D,
    const Population& Pop) {
    std::vector<Point> Matrix;

    // Берём dim ближайших соседей
    size_t dim = point.getDimension();
    std::vector<DistIndex> nearest = getKNearest(D, dim);

    Matrix.reserve(nearest.size());

    for (const auto& neighbor : nearest) {
        // Вектор от точки к соседу
        Point R = Pop[neighbor.getIndex()] - point;

        // Нормируем (превращаем в единичный вектор)
        long double len = R.norm();
        if (len > 1e-12) {
            for (size_t j = 0; j < dim; ++j) {
                R[j] /= len;
            }
        }

        Matrix.push_back(R);
    }

    return Matrix;
}

// Создание структуры соседей для всей популяции
std::vector<std::vector<DistIndex>> CreatePopulationSet(const Population& A) {
    size_t N = A.getSize();
    std::vector<std::vector<DistIndex>> DI(N);

    for (size_t i = 0; i < N; ++i) {
        DI[i].reserve(N - 1);
        for (size_t j = 0; j < N; ++j) {
            if (i == j) continue;
            long double dist = distance(A[i], A[j]);
            DI[i].emplace_back(j, dist);
        }
        std::sort(DI[i].begin(), DI[i].end(), CompareByDist);
    }
    return DI;
}

// Обновление структуры после изменения одной точки
void UpdatePopulation(Population& Pop,
    std::vector<std::vector<DistIndex>>& DI,
    const Point& point,
    size_t idx) {
    size_t N = Pop.getSize();

    // 0. Обновляем точку в популяции
    Pop[idx] = point;

    // Проверка размеров
    if (DI.empty() || DI[0].size() != N - 1) {
        std::cerr << "Error: Invalid DI structure" << std::endl;
        return;
    }

    // 1. Обновляем строку для изменённой точки idx
    DI[idx].clear();
    DI[idx].reserve(N - 1);
    for (size_t i = 0; i < N; ++i) {
        if (i == idx) continue;
        long double dist = distance(point, Pop[i]);
        DI[idx].emplace_back(i, dist);
    }
    std::sort(DI[idx].begin(), DI[idx].end(), CompareByDist);

    // 2. Обновляем строки для всех остальных точек
    for (size_t i = 0; i < N; ++i) {
        if (i == idx) continue;

        long double new_dist = distance(Pop[i], point);

        // Находим позицию старой записи об idx
        size_t old_pos = 0;
        bool found = false;
        for (size_t j = 0; j < DI[i].size(); ++j) {
            if (DI[i][j].getIndex() == idx) {
                old_pos = j;
                found = true;
                break;
            }
        }

        if (!found) {
            // Записи нет — добавляем новую
            DI[i].emplace_back(idx, new_dist);
            std::sort(DI[i].begin(), DI[i].end(), CompareByDist);
            continue;
        }

        // Удаляем старую запись
        DI[i].erase(DI[i].begin() + old_pos);

        // Вставляем новую запись в правильную позицию
        auto insert_pos = std::lower_bound(DI[i].begin(), DI[i].end(),
            DistIndex(idx, new_dist),
            CompareByDist);
        DI[i].insert(insert_pos, DistIndex(idx, new_dist));
    }
}

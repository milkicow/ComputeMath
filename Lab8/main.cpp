#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <chrono>

namespace mpi = boost::mpi;

class Transport {
    mpi::communicator *world_;
    double precision = 10e-5;
    double x_step_;
    double t_step_;

    double x_size_;
    double t_size_;
    std::vector<std::vector<double>>
        grid; // vector of x vectors for each moment t

    // std::vector<double> grid2_; // [x * t_size_ + t]

    void fill_x0();
    void fill_t0();
    double func(double x, double t);

  public:
    Transport(mpi::communicator *world, double h, double x_max, double t,
              double t_max)
        : world_(world), x_step_(h), t_step_(t),
        x_size_(x_max / h), t_size_(t_max / t) {

        grid.reserve(t_size_);
        for (auto &&vec : grid) {
            vec.reserve(x_size_);
        }
    }

    double u_x0(double x);
    double u_t0(double t);

    double corner2(double u, double u_next, double f);
    void grid_corner2();

    double corner(double u, double u_prev, double f);
    void grid_corner();

    double four_point(double u_prev, double u, double u_next, double f);
    void grid_four_point();

    double cross(double u_past, double u_prev, double u, double u_next,
                 double f);
    void grid_cross();

    void dump_grid();
};

void Transport::fill_x0() {
    std::vector<double> x0;
    double x = 0.0;

    for (auto &&it = 0; it < x_size_; ++it) {
        x0.emplace_back(u_x0(x));
        x += x_step_;
    }
    grid.emplace_back(std::move(x0));

    // std::vector<double> x0(x_size_);
    // grid.emplace_back(std::move(x0));
}

void Transport::fill_t0() {
    double t = 0;
    for (auto &&vec_it = 1; vec_it < t_size_; ++vec_it) {
        std::vector<double> vec_tmp(x_size_);
        vec_tmp[0] = u_t0(t);
        grid.push_back(std::move(vec_tmp));
        t += t_step_;
    }

    // double t = 0;
    // for (auto &&vec_it = 1; vec_it < t_size_; ++vec_it) {
    //     std::vector<double> vec_tmp(x_size_);
    //     vec_tmp[0] = 0;
    //     grid.push_back(std::move(vec_tmp));
    //     t += t_step_;
    // }
}

double Transport::func(double x, double t) { return t * x; }

void Transport::grid_corner() {
    fill_x0();
    fill_t0();

    for (int i = 0; i < grid.size() - 1; ++i) {
        for (int j = 0; j < grid[i].size() - 1; ++j) {
            grid[i + 1][j + 1] = corner(grid[i][j + 1], grid[i][j],
                                        func(x_step_ * j, t_step_ * i));
        }
    }
}

[[deprecated]] void Transport::grid_corner2() {
    fill_x0();
    fill_t0();

    for (int i = 0; i < grid.size() - 1; ++i) {
        for (int j = 0; j < grid[i].size() - 1; ++j) {
            grid[i + 1][j] = corner2(grid[i][j + 1], grid[i][j],
                                     func(x_step_ * j, t_step_ * i));
        }
        int j = grid[i].size() - 2;
        grid[i + 1][j + 1] =
            corner(grid[i][j + 1], grid[i][j], func(x_step_ * j, t_step_ * i));
    }
}

void Transport::grid_four_point() {
    fill_x0();
    fill_t0();

    for (int i = 0; i < grid.size() - 1; ++i) {
        for (int j = 1; j < grid[i].size() - 2; ++j) {
            grid[i + 1][j + 1] =
                four_point(grid[i][j], grid[i][j + 1], grid[i][j + 2],
                           func(x_step_ * j, t_step_ * i));
        }
        int j = grid[i].size() - 2;
        grid[i + 1][j + 1] =
            corner(grid[i][j + 1], grid[i][j], func(x_step_ * j, t_step_ * i));
    }
}

void Transport::grid_cross() {
    auto start = std::chrono::system_clock::now();
    fill_x0();
    fill_t0();

    assert(grid[0].size() % world_->size() == 0);
    auto size = grid[0].size() / world_->size();
    auto begin_j = world_->rank() * size;
    auto end_j = (world_->rank() + 1) * size;

    // [x * t_size_ + t] j == x, i == t
    // first layer
    for (int j = begin_j; j < end_j - 1; ++j) {
        grid[1][j + 1] =
            corner(grid[0][j + 1], grid[0][j], func(x_step_ * (j + 1), 0));
    }

    if (world_->rank() != world_->size() - 1) {
        world_->send(world_->rank() + 1, 0, grid[1][end_j - 1]);
    }
    double first_point_ = 0.0;
    //

    for (int i = 1; i < grid.size() - 1; ++i) {
        if (world_->rank() != 0) {
            world_->recv(world_->rank() - 1, 0, first_point_);
            grid[i + 1][begin_j] = first_point_;
        }

        for (int j = begin_j; j < end_j - 2; ++j) {
            grid[i + 1][j + 1] =
                cross(grid[i - 1][j + 1], grid[i][j], grid[i][j + 1],
                      grid[i][j + 2], func(x_step_ * (j + 1), t_step_ * i));
        }
        int j = end_j - 2;
        first_point_ = corner(grid[i][j + 1], grid[i][j],
                              func(x_step_ * (j + 1), t_step_ * i));
        grid[i + 1][j + 1] = first_point_;

        if (world_->rank() != world_->size() - 1) {
            world_->send(world_->rank() + 1, 0, first_point_);
        }
    }

    std::cout << world_->rank() << std::endl;
    auto end = std::chrono::system_clock::now();
    std::cout << world_->rank() <<" before bar time = " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << " secs" << std::endl;

    // world_->barrier();
    // Merge data
    if (world_->rank() != 0) {
        std::vector<double> send_data(grid.size() * size);
        for (int i = 0; i < grid.size(); ++i) {
            for (int j = 0; j < size; ++j) {
                send_data[(i * size) + j] = grid[i][j + begin_j];
            }
        }
        world_->send(0, 0, send_data);
    } else {
        for (int rank = 1; rank < world_->size(); ++rank) {
            std::vector<double> recv_data(grid.size() * size);
            world_->recv(rank, 0, recv_data);
            for (int i = 0; i < grid.size(); ++i) {
                for (int j = 0; j < size; ++j) {
                    grid[i][j + rank * size] = recv_data[(i * size) + j];
                }
            }
        }
    }
    world_->barrier();
}

void Transport::dump_grid() {
    std::cout << grid.size() << std::endl;

    std::ofstream dump_file{"output.txt"};
    dump_file << t_step_ << " " << x_step_ << std::endl;

    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[i].size(); ++j) {
            // std::cout << std::setw(5) << grid[i][j] << " ";
            dump_file << std::setw(5) << grid[i][j] << " ";
        }
        // std::cout << std::endl;
        dump_file << std::endl;
    }
}

// exaple from calc math: rect
double Transport::u_x0(double x) {
    return std::pow(x, 3) / 12;
}

// triangle
double Transport::u_t0(double t) {
    return std::pow(t, 3) / 12;
}

[[deprecated]] double Transport::corner2(double u, double u_next, double f) {
    return u * (1 + t_step_ / x_step_) - u_next * t_step_ / x_step_ +
           f * t_step_;
}

double Transport::corner(double u, double u_prev, double f) {
    return u * (1 - t_step_ / x_step_) + u_prev * t_step_ / x_step_ +
           f * t_step_;
}

double Transport::four_point(double u_prev, double u, double u_next, double f) {
    return u + f * t_step_ - (u_next - u_prev) * t_step_ / (2 * x_step_) +
           0.5 * t_step_ * t_step_ * (u_next - 2 * u + u_prev) /
               (x_step_ * x_step_);
}

double Transport::cross(double u_past, double u_prev, double u, double u_next,
                        double f) {
    return 2 * f * t_step_ - (u_next - u_prev) * t_step_ / x_step_ + u_past;
}

int main() {
    mpi::environment env;
    mpi::communicator world;

    Transport eq(&world, 10e-3, 10, 10e-3, 10);

    // eq.grid_four_point();
    // eq.grid_corner();
    auto start = std::chrono::system_clock::now();
    eq.grid_cross();
    auto end = std::chrono::system_clock::now();

    std::cout << world.rank() <<" time = " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << " secs" << std::endl;

    if (world.rank() == 0) {
        eq.dump_grid();
    }
    return 0;
}
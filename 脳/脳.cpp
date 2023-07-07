#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdarg>
#include <chrono>

#include "common.h"

template<typename T> T Sqr(const T &x) { return x * x; }

using namespace std;

struct Attendee {
    double x;
    double y;
    vector<double> tastes;
    FLD_BEGIN
        FLD(x) FLD(y) FLD(tastes)
    FLD_END
};

struct Problem {
    double room_width;
    double room_height;
    double stage_width;
    double stage_height;
    array<double, 2> stage_bottom_left;
    vector<int> musicians;
    vector<Attendee> attendees;
    FLD_BEGIN
        FLD(room_width) FLD(room_height) FLD(stage_width) FLD(stage_height) FLD(stage_bottom_left)
        FLD(musicians) FLD(attendees)
    FLD_END
};

struct Placement {
    double x, y;
    FLD_BEGIN FLD(x) FLD(y) FLD_END
};

struct Solution {
    vector<Placement> placements;
    FLD_BEGIN FLD(placements) FLD_END
};

// [musician][attendee]
vector<vector<int>> calc_visible(const Problem &p, const Solution &places) {
    // todo!!
    return {};
}

typedef double Weight;

vector<int> get_optimal_assignment(const vector<vector<Weight> > &a)
{
    const Weight inf = 1e30;
    int n = (int)a.size();
    if (n == 0) return vector<int>();
    int m = (int)a[0].size();
    if (n > m)
    {
        vector<vector<Weight> > at(m, vector<Weight>(n));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                at[j][i] = a[i][j];
        auto rest = get_optimal_assignment(at);
        vector<int> res(n, -1);
        for (int i = 0; i < m; i++) res[rest[i]] = i;
        return res;
    }
    vector<int> p(m + 1), way(m + 1);
    vector<Weight> u(n + 1), v(m + 1);
    vector<Weight> minv(m + 1);
    vector<char> used(m + 1);
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        fill(minv.begin(), minv.end(), inf);
        fill(used.begin(), used.end(), false);
        do
        {
            used[j0] = true;
            int i0 = p[j0], j1;
            Weight delta = inf;
            for (int j = 1; j <= m; ++j)
                if (!used[j])
                {
                    Weight cur = (i0 ? a[i0 - 1][j - 1] : 0) - u[i0] - v[j];
                    if (cur < minv[j])
                        minv[j] = cur, way[j] = j0;
                    if (minv[j] < delta)
                        delta = minv[j], j1 = j;
                }
            for (int j = 0; j <= m; ++j)
                if (used[j])
                    u[p[j]] += delta, v[j] -= delta;
                else
                    minv[j] -= delta;
            j0 = j1;
        } while (p[j0] != 0);
        do
        {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    vector<int> res(n, -1);
    for (int j = 1; j <= m; ++j)
        if (p[j] > 0)
            res[p[j] - 1] = j - 1;
    return res;
}

pair<Solution, double> solve_assignment(const Problem &p, const Solution &places) {
    auto visible = calc_visible(p, places);
    int n = (int)p.musicians.size();
    int m = (int)p.attendees.size();
    vector<vector<double>> mat(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int inst = p.musicians[i];
            double t = 0;
            for (int k = 0; k < m; k++) if (visible[j][k]) {
                double d2 = Sqr(places.placements[j].x - p.attendees[k].x) + Sqr(places.placements[j].y - p.attendees[k].y);
                t += ceil(1000000 * p.attendees[k].tastes[inst] / d2);
            }
            mat[i][j] = t;
        }
    auto ass = get_optimal_assignment(mat);
    Solution res;
    double score = 0;
    for (int i = 0; i < n; i++) {
        res.placements.push_back(places.placements[ass[i]]);
        score += mat[i][ass[i]];
    }
    return {res, score};
}

Solution get_some_placement(const Problem &p) {
    int n = (int)p.musicians.size();
    for (int i = 1; i <= n; n++) { // n rows
        int j = (n + i - 1) / i ; // n cols
        if (10 * i + 10 > p.stage_height || 10 * j + 10 > p.stage_width) continue;
        double dx = p.stage_width / (j + 1);
        double dy = p.stage_height / (i + 1);
        Solution res;
        for (int x = 0; x < j; x++)
            for (int y = 0; y < i; y++) {
                if ((int)res.placements.size() == n) break;
                res.placements.push_back({ dx * (x + 1), dy * (y + 1) });
            }
        return res;
    }
    return Solution();
}

void writeSolution(const Solution &sol, const char *tag, int problem_id) {
    auto f = fopen(format("../答/%d/%s.solution", problem_id, tag).c_str(), "wt");
    if (!f) exit(13);
    Json::FastWriter fw;
    fprintf(f, "%s", fw.write(serializeJson(sol)).c_str());
    fclose(f);
}

void solve(int problem_id) {
    Json::Value root, root_s;
    Problem p;
    auto fname = format("../問/%d.problem", problem_id);
    if (!readJsonFile(fname.c_str(), root)) {
        fprintf(stderr, "Invalid json 1!\n");
        exit(1);
    }
    if (!deserializeJson(p, root)) {
        fprintf(stderr, "Invalid json 3!\n");
        exit(1);
    }

    auto s0 = get_some_placement(p);
    if (s0.placements.empty()) {
        exit(2);
    }
    auto [s, score] = solve_assignment(p, s0);
    writeSolution(s, "test", problem_id);
    printf("%d musicians, %d attendees, %.0f x %.0f, score: %.3f", p.musicians.size(), p.attendees.size(), p.room_width, p.room_height, score);
}

int main() {
    solve(1);
    return 0;
}

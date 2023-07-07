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
    printf("%d musicians, %d attendees, %.0f x %.0f", p.musicians.size(), p.attendees.size(), p.room_width, p.room_height);
}

int main() {
    solve(1);
    return 0;
}

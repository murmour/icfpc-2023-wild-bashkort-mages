#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdarg>
#include <chrono>
#include <cstring>

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

struct T
{
	double x, y;
	T( double _x=0., double _y=0.) { x=_x; y=_y; }
};

bool is_valid( const Problem & problem, const Solution & sol )
{
	int n = (int)problem.musicians.size();
	if ((int)sol.placements.size() != n) return false;
	double sx = problem.stage_bottom_left[0];
	double sy = problem.stage_bottom_left[1];
	for (int i=0; i<n; i++)
	{
		double x = sol.placements[i].x;
		double y = sol.placements[i].y;
		if (!(sx <= x + 10 && x + 10 <= sx+problem.stage_width && sy <= y + 10 && y + 10 <= sy+problem.stage_height)) return false;
		for (int j=i+1; j<n; j++)
		{
			double dx = x - sol.placements[j].x;
			double dy = y - sol.placements[j].y;
			if (dx*dx + dy*dy < 100.0)
				return false;
		}
	}
	return true;
}

bool is_blocked( T A, T B, T C, double R )
{
	//return false;
	double dx1 = B.x - A.x, dy1 = B.y - A.y;
	double dx2 = C.x - A.x, dy2 = C.y - A.y;
	double cross1 = dx1*dx2 + dy1*dy2;
	if (cross1 <= 0.) return false;
	double dx3 = C.x - B.x, dy3 = C.y - B.y;
	double cross2 = -dx1*dx3 - dy1*dy3;
	if (cross2 <= 0.) return false;
	double S = abs( dx1*dy2 - dx2*dy1 );
	double d2 = dx1*dx1 + dy1*dy1;
	// R => h = S / sqrt(d2)
	// R*R => S*S / d2
	return R*R*d2 >= S*S;
}

vector< int > get_blocked_stupid( const Problem & problem, const Solution & sol, int mus_id )
{
	int n = (int)problem.musicians.size();
	int m = (int)problem.attendees.size();
	vector< int > res = vector< int >( m, 0 );

	T A = T( sol.placements[mus_id].x, sol.placements[mus_id].y );
	for (int i=0; i<m; i++)
	{
		T B = T( problem.attendees[i].x, problem.attendees[i].y );
		for (int j=0; j<n; j++)
			if (j != mus_id)
			{
				T C = T( sol.placements[j].x, sol.placements[j].y );
				if (is_blocked( A, B, C, 5. ))
				{
					res[i] = 1;
					break;
				}
			}
	}
	return res;
}

double pi = acos(-1.);

vector< int > get_blocked( const Problem & problem, const Solution & sol, int mus_id )
{
	vector< pair< double, int > > vec; // < angle, id >
	int n = (int)problem.musicians.size();
	T A = T( sol.placements[mus_id].x, sol.placements[mus_id].y );
	for (int i=0; i<n; i++)
		if (i!=mus_id)
		{
			T B = T( sol.placements[i].x, sol.placements[i].y );
			vec.push_back( make_pair( atan2( B.y-A.y, B.x-A.x ), i ) );
		}
	sort( vec.begin(), vec.end() );
	int m = (int)problem.attendees.size();
	vector< int > res = vector< int >( m, 0 );
	for (int i=0; i<m; i++)
	{
		T B = T( problem.attendees[i].x, problem.attendees[i].y );
		double angle = atan2( B.y-A.y, B.x-A.x );
		int tmp = lower_bound( vec.begin(), vec.end(), make_pair( angle, -1 ) ) - vec.begin();
		if (tmp==(int)vec.size()) tmp = 0;
		T C = T( sol.placements[vec[tmp].second].x, sol.placements[vec[tmp].second].y );
		if (is_blocked( A, B, C, 5. ))
		{
			res[i] = 1;
			continue;
		}
		int tmp2 = tmp+1;
		bool flag = false;
		while(true)
		{
			if (tmp2==(int)vec.size()) tmp2 = 0;
			if (tmp2==tmp) break;
			double delta = abs( vec[tmp2].first - angle );
			if (delta > pi) delta = abs( 2*pi - delta );
			if (delta*6 > pi) break;
			T C = T( sol.placements[vec[tmp2].second].x, sol.placements[vec[tmp2].second].y );
			if (is_blocked( A, B, C, 5. ))
			{
				res[i] = 1;
				flag = true;
				break;
			}
			tmp2++;
		}
		if (flag) continue;
		tmp2 = tmp-1;
		while(true)
		{
			if (tmp2==-1) tmp2 = (int)vec.size()-1;
			if (tmp2==tmp) break;
			double delta = abs( vec[tmp2].first - angle );
			if (delta > pi) delta = abs( 2*pi - delta );
			if (delta*6 > pi) break;
			T C = T( sol.placements[vec[tmp2].second].x, sol.placements[vec[tmp2].second].y );
			if (is_blocked( A, B, C, 5. ))
			{
				res[i] = 1;
				break;
			}
			tmp2--;
		}
	}
	return res;
}

vector< int > get_blocked2( const Problem & problem, const Solution & sol, int mus_id )
{
	int n = (int)problem.musicians.size();
	T A = T( sol.placements[mus_id].x, sol.placements[mus_id].y );
	vector< pair< double, int > > vecd; // < d2, id >

	for (int i=0; i<n; i++)
		if (i!=mus_id)
		{
			T B = T( sol.placements[i].x, sol.placements[i].y );
			double dx = A.x - B.x, dy = A.y - B.y;
			vecd.push_back( make_pair( dx*dx + dy*dy, i ) );
		}
	sort( vecd.begin(), vecd.end() );

	vector< vector< pair< double, int > > > vec; // < angle, id >
	double R = 9.;

	vector< pair< double, int > > vec_cur;
	vector< double > max_angle;
	for (int i=0; i<n; i++)
	{
		if (i==n-1 || vecd[i].first > 4*R*R)
		{
			if (vec_cur.size()>0)
			{
				sort( vec_cur.begin(), vec_cur.end() );
				vec.push_back( vec_cur );
				vec_cur.clear();
				max_angle.push_back( asin( 5./R ) + 0.0001 );
			}
			R *= 2;
		}
		if (i<n-1)
		{
			T B = T( sol.placements[vecd[i].second].x, sol.placements[vecd[i].second].y );
			vec_cur.push_back( make_pair( atan2( B.y-A.y, B.x-A.x ), vecd[i].second ) );
		}
	}

	int m = (int)problem.attendees.size();
	vector< int > res = vector< int >( m, 0 );
	for (int i=0; i<m; i++)
	{
		T B = T( problem.attendees[i].x, problem.attendees[i].y );
		double angle = atan2( B.y-A.y, B.x-A.x );
		for (int j=0; j<(int)vec.size(); j++)
		{
			int nn = vec[j].size();
			int tmp = lower_bound( vec[j].begin(), vec[j].end(), make_pair( angle, -1 ) ) - vec[j].begin();
			if (tmp==nn) tmp = 0;
			T C = T( sol.placements[vec[j][tmp].second].x, sol.placements[vec[j][tmp].second].y );
			if (is_blocked( A, B, C, 5. ))
			{
				res[i] = 1;
				break;
			}
			int tmp2 = tmp+1;
			bool flag = false;
			while(true)
			{
				if (tmp2==nn) tmp2 = 0;
				if (tmp2==tmp) break;
				double delta = abs( vec[j][tmp2].first - angle );
				if (delta > pi) delta = abs( 2*pi - delta );
				if (delta > max_angle[j]) break;
				T C = T( sol.placements[vec[j][tmp2].second].x, sol.placements[vec[j][tmp2].second].y );
				if (is_blocked( A, B, C, 5. ))
				{
					res[i] = 1;
					flag = true;
					break;
				}
				tmp2++;
			}
			if (flag) break;
			tmp2 = tmp-1;
			while(true)
			{
				if (tmp2==-1) tmp2 = nn-1;
				if (tmp2==tmp) break;
				double delta = abs( vec[j][tmp2].first - angle );
				if (delta > pi) delta = abs( 2*pi - delta );
				if (delta > max_angle[j]) break;
				T C = T( sol.placements[vec[j][tmp2].second].x, sol.placements[vec[j][tmp2].second].y );
				if (is_blocked( A, B, C, 5. ))
				{
					res[i] = 1;
					flag = true;
					break;
				}
				tmp2--;
			}
			if (flag) break;
		}
	}
	return res;
}

double get_score( const Problem & problem, const Solution & sol, bool use_ceil=true )
{
	int n = (int)problem.musicians.size();
	int m = (int)problem.attendees.size();
    assert(n == (int)sol.placements.size());
	double score = 0;

	for (int i=0; i<n; i++)
	{
		T A = T( sol.placements[i].x, sol.placements[i].y );
		vector< int > blocked = get_blocked2( problem, sol, i );
		for (int j=0; j<m; j++)
			if (blocked[j]==0)
			{
				T B = T( problem.attendees[j].x, problem.attendees[j].y );
				double dx = A.x - B.x, dy = A.y - B.y;
				double d2 = dx*dx + dy*dy;
				double tmp = 1'000'000 * problem.attendees[j].tastes[problem.musicians[i]];
                if (use_ceil)
				    score += ceil( tmp / d2 );
                else
                    score += tmp / d2;
			}
	}
	return score;
}

// [musician][attendee]
vector<vector<int>> calc_visible(const Problem &p, const Solution &places) {
	//cerr << "calc visible...";
	vector<vector<int>> res;
	int n = (int)p.musicians.size();
	int m = (int)p.attendees.size();
    for (int i=0; i<n; i++)
	{
		vector< int > tmp = get_blocked2( p, places, i );
		//vector< int > tmp2 = get_blocked_stupid( p, places, i );
		//if (tmp != tmp2) exit(666);
		for (int j=0; j<m; j++)
			tmp[j] = 1 - tmp[j];
		res.push_back( tmp );
	}
	//cerr << "ok\n";
	return res;
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
            int i0 = p[j0], j1 = 0;
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
	//cerr << "solve assignment\n";
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
            mat[i][j] = -t;
        }
    auto ass = get_optimal_assignment(mat);
    Solution res;
    double score = 0;
    for (int i = 0; i < n; i++) {
        res.placements.push_back(places.placements[ass[i]]);
        score += mat[i][ass[i]];
    }
    return {res, -score};
}

Solution get_regular_border_placement(const Problem & p, int mask = 15)
{
	int n = (int)p.musicians.size();
	double sx = p.stage_bottom_left[0];
	double sy = p.stage_bottom_left[1];
	Solution res;

	for (int side=0; side<4; side++)
		if ((mask>>side)&1)
		{
			if (side==0 || side==1)
			{
				for (int i=0; i*10+20 <= p.stage_width; i++)
				{
					double x = sx + 10 + i*10;
					double y = sy + (side==0 ? 10 : p.stage_height-10);
					bool flag = true;
					for (int j=0; j<(int)res.placements.size(); j++)
					{
						double dx = res.placements[j].x - x;
						double dy = res.placements[j].y - y;
						if (dx*dx + dy*dy < 100.)
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						res.placements.push_back( { x, y } );
						if ((int)res.placements.size() == n) return res;
					}
				}
			}

			if (side==2 || side==3)
			{
				for (int i=0; i*10+20 <= p.stage_height; i++)
				{
					double x = sx + (side==2 ? 10 : p.stage_width-10);
					double y = sy + 10 + i*10;
					bool flag = true;
					for (int j=0; j<(int)res.placements.size(); j++)
					{
						double dx = res.placements[j].x - x;
						double dy = res.placements[j].y - y;
						if (dx*dx + dy*dy < 100.)
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						res.placements.push_back( { x, y } );
						if ((int)res.placements.size() == n) return res;
					}
				}
			}
		}
	while ((int)res.placements.size() < n)
	{
		double x = sx + (p.stage_width-20.) * rand()/(RAND_MAX-1) + 10.;
		double y = sy + (p.stage_height-20.) * rand()/(RAND_MAX-1) + 10.;
		bool flag = true;
		for (int j=0; j<(int)res.placements.size(); j++)
		{
			double dx = res.placements[j].x - x;
			double dy = res.placements[j].y - y;
			if (dx*dx + dy*dy < 100.)
			{
				flag = false;
				break;
			}
		}
		if (flag) res.placements.push_back( { x, y } );
	}
	return res;
}

Solution get_border_placement(const Problem & p, int mask = 15)
{
	int n = (int)p.musicians.size();
	double sx = p.stage_bottom_left[0];
	double sy = p.stage_bottom_left[1];
	Solution res;
	for (int i=0; i<n; i++)
	{
		int iters = 0;
		while(true)
		{
			iters++;
			if (iters > 1000)
			{
				double x = sx + (p.stage_width-20.) * rand()/(RAND_MAX-1) + 10.;
				double y = sy + (p.stage_height-20.) * rand()/(RAND_MAX-1) + 10.;
				bool flag = true;
				for (int j=0; j<i; j++)
				{
					double dx = res.placements[j].x - x;
					double dy = res.placements[j].y - y;
					if (dx*dx + dy*dy < 100.)
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					res.placements.push_back( { x, y } );
					break;
				}
				else continue;
			}
			int side = rand()%4;
			if ( ((mask >> side)&1)==0 ) continue;
			double x, y;
			if (side==0 || side==1)
				x = sx + (p.stage_width-20.) * rand()/(RAND_MAX-1) + 10.;
			else if (side==2)
				x = sx + 10.;
			else x = sx + p.stage_width - 10.;
			if (side==2 || side==3)
				y = sy + (p.stage_height-20.) * rand()/(RAND_MAX-1) + 10.;
			else if (side==0)
				y = sy + 10.;
			else y = sy + p.stage_height - 10.;
			bool flag = true;
			for (int j=0; j<i; j++)
			{
				double dx = res.placements[j].x - x;
				double dy = res.placements[j].y - y;
				if (dx*dx + dy*dy < 100.)
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				res.placements.push_back( { x, y } );
				break;
			}
		}
	}
	//cerr << "got border placement\n";
	return res;
}

Solution get_two_row_border_placement(const Problem & p, int mask = 15)
{
	int n = (int)p.musicians.size();
	double sx = p.stage_bottom_left[0];
	double sy = p.stage_bottom_left[1];
	Solution res;
	for (int i=0; i<n; i++)
	{
		int iters = 0;
		while(true)
		{
			iters++;
			if (iters > 1000)
			{
				double x = sx + (p.stage_width-20.) * rand()/(RAND_MAX-1) + 10.;
				double y = sy + (p.stage_height-20.) * rand()/(RAND_MAX-1) + 10.;
				bool flag = true;
				for (int j=0; j<i; j++)
				{
					double dx = res.placements[j].x - x;
					double dy = res.placements[j].y - y;
					if (dx*dx + dy*dy < 100.)
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					res.placements.push_back( { x, y } );
					break;
				}
				else continue;
			}
			int side = rand()%4;
			if ( ((mask >> side)&1)==0 ) continue;
			double x, y;
			if (side==0 || side==1)
				x = sx + (p.stage_width-20.) * rand()/(RAND_MAX-1) + 10.;
			else if (side==2)
				x = sx + 10. + 10.*(rand()&1);
			else x = sx + p.stage_width - 10. - 10.*(rand()&1);
			if (side==2 || side==3)
				y = sy + (p.stage_height-20.) * rand()/(RAND_MAX-1) + 10.;
			else if (side==0)
				y = sy + 10. + 10.*(rand()&1);
			else y = sy + p.stage_height - 10. - 10.*(rand()&1);
			bool flag = true;
			for (int j=0; j<i; j++)
			{
				double dx = res.placements[j].x - x;
				double dy = res.placements[j].y - y;
				if (dx*dx + dy*dy < 100.)
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				res.placements.push_back( { x, y } );
				break;
			}
		}
	}
	//cerr << "got border placement\n";
	return res;
}

Solution get_random_placement(const Problem & p)
{
	//cerr << "get random placement\n";
	int n = (int)p.musicians.size();
	double sx = p.stage_bottom_left[0];
	double sy = p.stage_bottom_left[1];
	Solution res;
	for (int i=0; i<n; i++)
	{
		while(true)
		{
			double x = sx + (p.stage_width-20.) * rand()/(RAND_MAX-1) + 10.;
			double y = sy + (p.stage_height-20.) * rand()/(RAND_MAX-1) + 10.;
			bool flag = true;
			for (int j=0; j<i; j++)
			{
				double dx = res.placements[j].x - x;
				double dy = res.placements[j].y - y;
				if (dx*dx + dy*dy < 100.)
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				res.placements.push_back( { x, y } );
				break;
			}
		}
	}
	return res;
}

Solution get_some_placement(const Problem &p) {
	cerr << "get some placement\n";
    int n = (int)p.musicians.size();
	double sx = p.stage_bottom_left[0];
	double sy = p.stage_bottom_left[1];
    for (int i = 1; i <= n; i++) { // n rows
        int j = (n + i - 1) / i ; // n cols
        if (10 * i + 10 > p.stage_height || 10 * j + 10 > p.stage_width) continue;
        double dx = p.stage_width / (j + 1);
        double dy = p.stage_height / (i + 1);
        Solution res;
        for (int x = 0; x < j; x++)
            for (int y = 0; y < i; y++) {
                if ((int)res.placements.size() == n) break;
                res.placements.push_back({ sx + dx * (x + 1), sy + dy * (y + 1) });
            }
        return res;
    }
    return Solution();
}

void writeSolution(const Solution &sol, string fname) {
    //auto f = fopen(format("../答/%d/%s.solution", problem_id, tag).c_str(), "wt");
    if (!fname.empty()) {
        auto f = fopen(fname.c_str(), "wt");
        if (!f) exit(13);
        Json::FastWriter fw;
        fprintf(f, "%s", fw.write(serializeJson(sol)).c_str());
        fclose(f);
    } else {
        Json::FastWriter fw;
        printf("%s", fw.write(serializeJson(sol)).c_str());
    }
}

void solve(const string &infile, int timeout, const string &solver, const string &fname) {
    Json::Value root, root_s;
    Problem p;
    if (!readJsonFile(infile.c_str(), root)) {
        fprintf(stderr, "Invalid json 1!\n");
        exit(1);
    }
    if (!deserializeJson(p, root)) {
        fprintf(stderr, "Invalid json 3!\n");
        exit(1);
    }

	fprintf(stderr, "input %s, solver %s, %d musicians, %d attendees, %.0lf x %.0lf\n", infile.c_str(), solver.c_str(), (int)p.musicians.size(), (int)p.attendees.size(), p.room_width, p.room_height );

	double best_score = 0.;
	int iters = 0;
	int start_time = clock();
    Solution best_solution;
	while(true)
	{
		int cur_time = clock();
		if (cur_time - start_time > timeout*CLOCKS_PER_SEC) break;
        Solution s0;
        if (solver == "two_row")
		    s0 = get_two_row_border_placement(p, rand()%16);
        else if (solver == "border")
            s0 = get_border_placement(p, rand() % 16);
        else if (solver == "regular")
            s0 = get_regular_border_placement(p, rand() % 16);
        else {
            fprintf(stderr, "Invalid solver: %s\n", solver.c_str());
			exit(10);
        }

		if (s0.placements.empty()) {
			exit(2);
		}
		if (!is_valid(p,s0))
		{
			fprintf(stderr, "Invalid placement!\n");
			exit(3);
		}
		auto [s, score] = solve_assignment(p, s0);
		iters++;
		if (score > best_score)
		{
			best_score = score;
            best_solution = s;
			fprintf(stderr, "iters: %d score: %.3lf\n", iters, score);
            if (!fname.empty()) writeSolution(s, fname);
			//double my_score = get_score(p,s);
			//printf("my score %.3lf\n", my_score);
		}
	}
    if (best_solution.placements.empty()) {
        exit(11);
    }
    writeSolution(best_solution, fname);
}

struct ArgParser {
	int argc;
	char **argv;
	char * get_arg(const char * name) {
		for (int i = 1; i + 1 < argc; i += 2)
			if (strcmp(argv[i], name) == 0)
				return argv[i + 1];
		return nullptr;
	}
};

int main(int argc, char *argv[]) {
    ArgParser args = { argc, argv };

    string in_file = "";
    if (auto p = args.get_arg("-pp"))
    {
        in_file = p;
    } else {
        if (auto p = args.get_arg("-p"))
            in_file = format("../問/%s.problem", p);
    }

    int timeout = 120;
    if (auto p = args.get_arg("-timeout"))
        sscanf(p, "%d", &timeout);

    string solver = "regular";
	if (auto s = args.get_arg("-s")) {
        solver = s;
    }

    string fname = "";
	if (auto s = args.get_arg("-out")) {
        fname = s;
    }

    if (auto p = args.get_arg("-score")) {
        Json::Value root;
        Solution sol;
        if (!readJsonFile(p, root)) { return 1; }
        if (!deserializeJson(sol, root)) { return 2; }
        Problem prob;
        if (!readJsonFile(in_file.c_str(), root)) {
            fprintf(stderr, "Invalid json 1!\n");
            exit(1);
        }
        if (!deserializeJson(prob, root)) {
            fprintf(stderr, "Invalid json 3!\n");
            exit(1);
        }
		if (!is_valid(prob, sol)) {
			fprintf(stderr, "Invalid solution!\n");
            exit(10);
		}
        double score = get_score(prob, sol);
        printf("%.3lf", score);
        return 0;
    }

    solve(in_file, timeout, solver, fname);
	//for (int i=1; i<=45; i++)
	//	solve(i);
    return 0;
}

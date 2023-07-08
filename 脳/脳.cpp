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
#include <regex>

#include "common.h"
#include "geom2d.h"

template<typename T> T Sqr(const T &x) { return x * x; }

using namespace std;

int problem_id = -1;
bool new_scoring = false;

struct Attendee {
    double x;
    double y;
    vector<double> tastes;
    FLD_BEGIN
        FLD(x) FLD(y) FLD(tastes)
    FLD_END
};

struct Pillar {
	array<double, 2> center;
	double radius;
	FLD_BEGIN FLD(center) FLD(radius) FLD_END
};

struct Problem {
    double room_width;
    double room_height;
    double stage_width;
    double stage_height;
    array<double, 2> stage_bottom_left;
    vector<int> musicians;
    vector<Attendee> attendees;
	vector<Pillar> pillars;
    FLD_BEGIN
        FLD(room_width) FLD(room_height) FLD(stage_width) FLD(stage_height) FLD(stage_bottom_left)
        FLD(musicians) FLD(attendees) FLD(pillars)
    FLD_END
};

struct Placement {
    double x, y;
    FLD_BEGIN FLD(x) FLD(y) FLD_END
};

struct Solution {
    vector<Placement> placements;
	double score = -1;
    FLD_BEGIN FLD(placements) FLD(score, -1) FLD_END
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
		if (!(sx <= x + 10 && x + 10 <= sx+problem.stage_width && sy <= y + 10 && y + 10 <= sy+problem.stage_height))
		{
			fprintf(stderr, "mus %d is out of stage: (%.3f, %.3f)\n", i, x, y);
			return false;
		}
		for (int j=i+1; j<n; j++)
		{
			double dx = x - sol.placements[j].x;
			double dy = y - sol.placements[j].y;
			if (dx*dx + dy*dy < 100.0) {
				double dist = hypot(dx, dy);
				fprintf(stderr, "mus %d and %d are too close: %.8f\n", i, j, dist);
				return false;
			}
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

	auto q = vector<double>(n);
	if (new_scoring) {
		for (int i = 0; i < n; i++) {
			double t = 1;
			for (int j = 0; j < n; j++) if (i != j && problem.musicians[i] == problem.musicians[j]) t += 1.0 / hypot(sol.placements[i].x - sol.placements[j].x, sol.placements[i].y - sol.placements[j].y);
			q[i] = t;
		}
	} else {
		for (int i = 0; i < n; i++) q[i] = 1.0;
	}

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
				    score += ceil(q[i] * ceil( tmp / d2 ));
                else
                    score += q[i] * tmp / d2;
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

Solution solve_assignment(const Problem &p, const Solution &places) {
	//cerr << "solve assignment\n";
    auto visible = calc_visible(p, places);
    int n = (int)p.musicians.size();
    int m = (int)p.attendees.size();
	/*
	auto q = vector<double>(n);
	if (new_scoring) {
		for (int i = 0; i < n; i++) {
			double t = 1;
			for (int j = 0; j < n; j++) if (i != j && p.musicians[i] == p.musicians[j]) t += 1.0 / hypot(places.placements[i].x - places.placements[j].x, places.placements[i].y - places.placements[j].y);
			q[i] = t;
		}
	} else {
		for (int i = 0; i < n; i++) q[i] = 1.0;
	}
	*/
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
	if (new_scoring)
		res.score = get_score(p, res);
	else
		res.score = -score;
    return res;
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

Solution get_compact_placement(const Problem &p, int xmode = 0, int ymode = 0) {
	int n = (int)p.musicians.size();
	int k = (int)ceil(sqrt(n));
	// k * k
	Solution res;
	double a = 10;
	double b = p.stage_width - k * 10;
	double x0 = a + (b - a) / 2 * xmode;
	b = p.stage_height - k * 10;
	double y0 = a + (b - a) / 2 * ymode;
	for (int i = 0; i < k; i++)
		for (int j = 0; j < k; j++) {
			if ((int)res.placements.size() >= n) break;
			double x = p.stage_bottom_left[0] + x0 + i * 10;
			double y = p.stage_bottom_left[1] + y0 + j * 10;
			res.placements.push_back({x, y});
		}
	return res;
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

inline double score_f_no_ceil(double x1, double y1, double x2, double y2, double taste) {
	double d2 = Sqr(x1 - x2) + Sqr(y1 - y2);
	double tmp = 1'000'000 * taste;
	return tmp / d2;
}

Solution wiggle(const Problem &p, const Solution &sol) {
    auto visible = calc_visible(p, sol);
	int n = (int)p.musicians.size();
	int m = (int)p.attendees.size();
	const double eps = 1e-9;

	double sx1 = p.stage_bottom_left[0] + 10;
	double sx2 = sx1 + p.stage_width - 20;
	double sy1 = p.stage_bottom_left[1] + 10;
	double sy2 = sy1 + p.stage_height - 20;

	Solution res = sol;
	auto &cur = res.placements;
	double dstep = 0.01;
	double step = 10;

	for (int big_iter = 0; big_iter < 10; big_iter++) {
		bool changed = false;
		for (int iter = 0; iter < 500; iter++) {
			double best_delta = 0;
			int best_idx = -1;
			//double best_dist;
			Placement best_pos;

			for (int i = 0; i < n; i++) {
				double x0 = cur[i].x;
				double y0 = cur[i].y;
				double base = 0;
				double dx = 0;
				double dy = 0;
				for (int j = 0; j < m; j++) if (visible[i][j]) {
					double taste = p.attendees[j].tastes[p.musicians[i]];
					base += score_f_no_ceil(x0, y0, p.attendees[j].x, p.attendees[j].y, taste);
					dx += score_f_no_ceil(x0 + dstep, y0, p.attendees[j].x, p.attendees[j].y, taste);
					dy += score_f_no_ceil(x0, y0 + dstep, p.attendees[j].x, p.attendees[j].y, taste);
				}
				Point grad(dx - base, dy - base);
				if (x0 <= sx1 + eps) grad.x = max(grad.x, 0.0);
				if (y0 <= sy1 + eps) grad.y = max(grad.y, 0.0);
				if (x0 >= sx2 - eps) grad.x = min(grad.x, 0.0);
				if (y0 >= sy2 - eps) grad.y = min(grad.y, 0.0);
				if (grad.len2() == 0) continue; // probably noone is visible
				auto dir = grad.normalized();
				// get max distance that we can move
				double dist = step;
				if (dir.x > eps)
					dist = min(dist, (sx2 - x0) / dir.x);
				if (dir.x < -eps)
					dist = min(dist, (sx1 - x0) / dir.x);
				if (dir.y > eps)
					dist = min(dist, (sy2 - y0) / dir.y);
				if (dir.y < -eps)
					dist = min(dist, (sy1 - y0) / dir.y);
				Point p0(x0, y0);
				auto line = Line::PP(p0, p0 + dir);
				for (int j = 0; j < n; j++) if (i != j) {
					auto h = abs(line.at({ cur[j].x, cur[j].y }));
					if (h >= 10) continue;
					auto v = p0.to({ cur[j].x, cur[j].y });
					double t = v.dot(dir);
					if (t < 0) continue;
					double d = t - sqrt(100 - h * h);
					dist = min(dist, d);
				}
				dist = max(0.0, dist - 1e-6);
				auto tgt = p0 + dir * dist;
				double new_score = 0;
				for (int j = 0; j < m; j++) if (visible[i][j]) {
					double taste = p.attendees[j].tastes[p.musicians[i]];
					new_score += score_f_no_ceil(tgt.x, tgt.y, p.attendees[j].x, p.attendees[j].y, taste);
				}
				double delta = new_score - base;
				if (delta > best_delta) {
					best_delta = delta;
					best_idx = i;
					best_pos = {tgt.x, tgt.y};
					//best_dist = dist;
				}
			}
			if (best_idx != -1) {
				changed = true;
				//fprintf(stderr, "delta = %.3f, idx = %d, dist = %.3f, tgt = (%.3f, %.3f)\n", best_delta, best_idx, best_dist, best_pos.x, best_pos.y);
				cur[best_idx] = best_pos;
			} else {
				break;
			}
			//if (best_dist < 0.1 * step) break;
		}
		if (changed) visible = calc_visible(p, res);
		step = step * 0.5;
	}
	return res;
}

inline double score_f(double x1, double y1, double x2, double y2, double taste) {
	double d2 = Sqr(x1 - x2) + Sqr(y1 - y2);
	double tmp = 1'000'000 * taste;
	return ceil(tmp / d2);
}

void print_stats(const Problem &p) {
	double cx = p.stage_bottom_left[0] + 0.5 * p.stage_width;
	double cy = p.stage_bottom_left[1] + 0.5 * p.stage_height;
	int n_instruments = (int)p.attendees[0].tastes.size();
	int m = (int)p.attendees.size();
	for (int i = 0; i < n_instruments; i++) {
		double score = 0;
		for (int j = 0; j < m; j++)
			score += score_f(cx, cy, p.attendees[j].x, p.attendees[j].y, p.attendees[j].tastes[i]);
		printf("%d: %.0f\n", i, score);
	}
}

vector<double> get_instrument_scores(const Problem &p) {
	double cx = p.stage_bottom_left[0] + 0.5 * p.stage_width;
	double cy = p.stage_bottom_left[1] + 0.5 * p.stage_height;
	int n_instruments = (int)p.attendees[0].tastes.size();
	int m = (int)p.attendees.size();
	vector<double> res;
	for (int i = 0; i < n_instruments; i++) {
		double score = 0;
		for (int j = 0; j < m; j++)
			score += score_f(cx, cy, p.attendees[j].x, p.attendees[j].y, p.attendees[j].tastes[i]);
		res.push_back(score);
	}
	return res;
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
		else if (solver == "wiggle") {
			//if (iters > 0) break;
			//s0 = get_random_placement(p);
			//s0 = get_border_placement(p, iters % 16);
			s0 = get_regular_border_placement(p, iters % 16);
			s0 = solve_assignment(p, s0);
			for (int i = 0; i < 1; i++) {
				//auto before = s0.score;
				s0 = wiggle(p, s0);
				s0 = solve_assignment(p, s0);
				//auto after = s0.score;
				//fprintf(stderr, "wiggle = %.0f\n", after - before);
			}
		}
		else if (solver == "compact") {
			if (iters >= 9) break;
			s0 = get_compact_placement(p, iters % 3, iters / 3);
		}
		else if (solver == "stats") {
			print_stats(p);
			break;
		}
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
		auto s = solve_assignment(p, s0);
		iters++;
		if (s.score > best_score)
		{
			best_score = s.score;
            best_solution = s;
			fprintf(stderr, "p%d: iters: %d score: %.3lf\n", problem_id, iters, s.score);
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

	regex prob_id("/(\\d+)\\.problem");
	smatch match;
	if (regex_search(in_file, match, prob_id)) {
		string s = match[1].str();
		sscanf(s.c_str(), "%d", &problem_id);
		if (problem_id >= 56) new_scoring = true;
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
	if (auto s = args.get_arg("-tag")) {
		fname = format("../答/%d/%s.solution", problem_id, s);
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

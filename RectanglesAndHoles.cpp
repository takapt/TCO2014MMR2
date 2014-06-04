// #define NDEBUG

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <map>
#include <utility>
#include <set>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <complex>
#include <stack>
#include <queue>
#include <numeric>
#include <list>
#include <iomanip>
#include <fstream>
#include <bitset>

using namespace std;

#define foreach(it, c) for (__typeof__((c).begin()) it=(c).begin(); it != (c).end(); ++it)
template <typename T> void print_container(ostream& os, const T& c) { const char* _s = " "; if (!c.empty()) { __typeof__(c.begin()) last = --c.end(); foreach (it, c) { os << *it; if (it != last) os << _s; } } }
template <typename T> ostream& operator<<(ostream& os, const vector<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const set<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const multiset<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const deque<T>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const map<T, U>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const pair<T, U>& p) { os << "(" << p.first << ", " << p.second << ")"; return os; }

template <typename T> void print(T a, int n, const string& split = " ") { for (int i = 0; i < n; i++) { cout << a[i]; if (i + 1 != n) cout << split; } cout << endl; }
template <typename T> void print2d(T a, int w, int h, int width = -1, int br = 0) { for (int i = 0; i < h; ++i) { for (int j = 0; j < w; ++j) { if (width != -1) cout.width(width); cout << a[i][j] << ' '; } cout << endl; } while (br--) cout << endl; }
template <typename T> void input(T& a, int n) { for (int i = 0; i < n; ++i) cin >> a[i]; }
#define dump(v) (cerr << #v << ": " << v << endl)

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define clr(a, x) memset(a, x, sizeof(a))
#define sz(a) ((int)(a).size())
#define mp(a, b) make_pair(a, b)
#define ten(n) ((long long)(1e##n))

template <typename T, typename U> void upmin(T& a, const U& b) { a = min<T>(a, b); }
template <typename T, typename U> void upmax(T& a, const U& b) { a = max<T>(a, b); }
template <typename T> void uniq(T& a) { sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end()); }
template <class T> string to_s(const T& a) { ostringstream os; os << a; return os.str(); }
template <class T> T to_T(const string& s) { istringstream is(s); T res; is >> res; return res; }
void fast_io() { cin.tie(0); ios::sync_with_stdio(false); }
bool in_rect(int x, int y, int w, int h) { return 0 <= x && x < w && 0 <= y && y < h; }


typedef pair<int, int> pint;
typedef long long ll;

const int DX[] = { 0, 1, 0, -1 };
const int DY[] = { 1, 0, -1, 0 };


#ifdef _MSC_VER
#include <Windows.h>
#else
#include <sys/time.h>
#endif
class Timer
{
    typedef double time_type;
    typedef unsigned int skip_type;

private:
    time_type start_time;
    time_type elapsed;

#ifdef _MSC_VER
    time_type get_ms() { return (time_type)GetTickCount64() / 1000; }
#else
    time_type get_ms() { struct timeval t; gettimeofday(&t, NULL); return (time_type)t.tv_sec * 1000 + (time_type)t.tv_usec / 1000; }
#endif

public:
    Timer() {}

    void start() { start_time = get_ms(); }
    time_type get_elapsed() { return elapsed = get_ms() - start_time; }
};


enum Dir
{
    UP,
    RIGHT,
    DOWN,
    LEFT,

    NA,
};
const string dir_s[] = { "UP", "RIGHT", "DOWN", "LEFT", "NA" };
Dir to_dir(const string& s)
{
    int i = find(dir_s, dir_s + 4, s) - dir_s;
    assert(0 <= i && i < 4);
    return Dir(i);
}
Dir rev_dir(Dir dir)
{
    return Dir((dir + 2) % 4);
}
struct Pos
{
    int x, y;
    Pos(int x, int y)
        : x(x), y(y)
    {
    }
    Pos()
        : x(0), y(0)
    {
    }

    bool operator==(const Pos& other) const
    {
        return x == other.x && y == other.y;
    }
    bool operator !=(const Pos& other) const
    {
        return x != other.x || y != other.y;
    }

    void operator+=(const Pos& other)
    {
        x += other.x;
        y += other.y;
    }
    void operator-=(const Pos& other)
    {
        x -= other.x;
        y -= other.y;
    }

    Pos operator+(const Pos& other) const
    {
        Pos res = *this;
        res += other;
        return res;
    }
    Pos operator-(const Pos& other) const
    {
        Pos res = *this;
        res -= other;
        return res;
    }
    Pos operator*(int a) const
    {
        return Pos(x * a, y * a);
    }

    bool operator<(const Pos& other) const
    {
        if (x != other.x)
            return x < other.x;
        else
            return y < other.y;
    }
};
Pos operator*(int a, const Pos& cell)
{
    return cell * a;
}
Pos to_cell(Dir dir)
{
    assert(0 <= dir && dir < 4);
    return Pos(DX[dir], DY[dir]);
}
namespace std
{
    ostream& operator<<(ostream& os, Dir dir)
    {
        os << dir_s[dir];
        return os;
    }

    ostream& operator<<(ostream& os, const Pos& c)
    {
        char buf[256];
        sprintf(buf, "(%d, %d)", c.x, c.y);
        os << buf;
        return os;
    }
}


const int MAX_RANGE = ten(6);
const int LOWER = -MAX_RANGE;
const int UPPER = MAX_RANGE;
const int MAX_SIDE_LEN = 1000;

bool check_pos(const Pos& pos)
{
    return max(abs(pos.x), abs(pos.y)) <= MAX_RANGE;
}

class Rect
{
public:
    Rect(int a, int b, int index)
        : w(a), h(b), index_(index), rotated_(false), pos_(-1919810, -114514)
    {
    }
    Rect()
        : index_(-810114514){}

    int width() const { return w; }
    int height() const { return h; }
    int index() const { return index_; }
    bool rotated() const { return rotated_; }
    void rotate()
    {
        swap(w, h);
        rotated_ ^= true;
    }

    int long_len() const { return max(w, h); }
    int short_len() const { return min(w, h); }

    void width_to_long()
    {
        if (w < h)
            rotate();
    }
    void height_to_long()
    {
        if (w > h)
            rotate();
    }

    const Pos& pos() const { return pos_; }
    void set_pos(const Pos& pos)
    {
        pos_ = pos;
    }

private:
    int w, h;
    int index_;
    bool rotated_;
    Pos pos_;
};


class WeenCreater;
class Ween
{
public:
    Ween()
        : pos_(Pos(0, 0))
    {
    }

    const Pos& pos() const { return pos_; }
    void set_pos(const Pos& pos)
    {
        Pos move = pos - pos_;
        pos_ = pos;

        for (auto& rect : rects_)
            rect.set_pos(rect.pos() + move);
    }

    vector<Rect> rects() const
    {
        vector<Rect> res;
        for (auto& rect : rects_)
            res.push_back(rect);

        for (auto& rect : res)
            rect.set_pos(pos_ + rect.pos());
        return res;
    }

private:
    // positions are relative to pos_
    vector<Rect> rects_;

    Pos pos_; // absolute pos in solution plane

    friend WeenCreater;
};
class WeenCreater
{
public:
    WeenCreater(){}
    WeenCreater(const vector<Rect>& rects_)
        : rects(rects_)
    {
        vector<int> sorted;
        rep(i, rects.size())
            sorted.push_back(i);
        sort(all(sorted), [&](int i, int j) { return rects[i].long_len() < rects[j].long_len(); });
        sorted.erase(sorted.begin() + rects.size() * 3 / 4, sorted.end());

        for (int i : sorted)
            outer_cand.push_back(sorted[i]);
        sort(all(outer_cand), [&](int i, int j) { return rects[i].long_len() < rects[j].long_len(); });

        // TODO: (long, index), (short, index)の両方をinner_candに突っ込んで決めたほうがいいか？
//         rep(i, rects.size())
        for (int i : sorted)
        {
            if (double(rects[i].long_len()) / rects[i].short_len() > 2)
                inner_cand.insert(pint(rects[i].long_len(), i));
        }
    }

    bool create(Ween& ween, vector<bool>& used)
    {
        ween = Ween();

        Rect outer[4];
        //
        rep(i, 4)
        {
            while (!outer_cand.empty() && used[rects[outer_cand.back()].index()])
                outer_cand.pop_back();
            if (outer_cand.empty())
                return false;

            outer[i] = rects[outer_cand.back()];
            outer_cand.pop_back();
        }
        rep(i, 4)
            used[outer[i].index()] = true;

        // 
        Rect& bottom = outer[0];
        bottom.width_to_long();

        Rect& top = outer[1];
        top.width_to_long();

        Rect& left = outer[2];
        left.height_to_long();

        Rect& right = outer[3];
        right.height_to_long();

        bottom.set_pos(Pos(0, -bottom.height()));
        top.set_pos(Pos(0, right.height()));
        left.set_pos(Pos(-left.width(), 0));
        ween.rects_.push_back(bottom);
        ween.rects_.push_back(top);
        ween.rects_.push_back(left);
        create_recur(right, top.width() + right.width() + 1, right.height(), MoveDir::X, ween, used);

#ifndef NDEBUG
        rep(i, 4)
        {
            assert(used[ween.rects_[i].index()]);
            assert(check_pos(ween.rects_[i].pos()));
        }
#endif

        return true;
    }
private:
    enum MoveDir { X, Y };
    void create_recur(Rect ween_rect, int width, int height, MoveDir move_dir, Ween& ween, vector<bool>& used)
    {
        if (move_dir == X)
        {
            const int next_width_range = width - ween_rect.width() - 1;
            assert(next_width_range > 0);
            Rect next_ween_rect;
            if (!select_ween_rect(next_width_range, height - 2, next_ween_rect, used))
            {
                ween_rect.set_pos(Pos(1, 0));
                used[ween_rect.index()] = true;
                ween.rects_.push_back(ween_rect);
                return;
            }

            next_ween_rect.width_to_long();
            ween_rect.set_pos(Pos(next_ween_rect.width(), 0));
            used[ween_rect.index()] = true;
            ween.rects_.push_back(ween_rect);

            create_recur(next_ween_rect, next_ween_rect.width(), height, MoveDir::Y, ween, used);
        }
        else if (move_dir == Y)
        {
            const int next_height_range = height - ween_rect.height() - 1;
            assert(next_height_range > 0);
            Rect next_ween_rect;
            if (!select_ween_rect(next_height_range, width - 2, next_ween_rect, used))
            {
                ween_rect.set_pos(Pos(0, 1));
                used[ween_rect.index()] = true;
                ween.rects_.push_back(ween_rect);
                return;
            }

            next_ween_rect.height_to_long();
            ween_rect.set_pos(Pos(0, next_ween_rect.height()));
            used[ween_rect.index()] = true;
            ween.rects_.push_back(ween_rect);

            create_recur(next_ween_rect, width, next_ween_rect.height(), MoveDir::X, ween, used);
        }
        else
            abort();
    }

    bool select_ween_rect(int long_range, int short_range, Rect& selected_rect, vector<bool>& used)
    {
        for (auto it = inner_cand.lower_bound(pint(long_range, -1)); it != inner_cand.end(); )
        {
            auto next = it;
            ++next;

            auto& rect = rects[it->second];
            assert(rect.long_len() <= long_range);

            if (used[rect.index()])
                inner_cand.erase(it);
            else if (rect.short_len() <= short_range)
            {
                selected_rect = rect;
                return true;
            }

            it = next;
        }
        return false;
    }

private:
    vector<Rect> rects;

    vector<int> outer_cand;
    set<pint, greater<pint>> inner_cand; // pint(long_len, rects index)
};

class Score
{
public:
    Score(int holes, ll area)
        : holes_(holes), area_(area)
    {
    }

    int holes() const { return holes_; }
    ll area() const { return area_; }
    ll score() const { return area_ * holes_ * holes_; }

    ll darea;
    ll dscore;
    double hole_cost;
    int exp_holes;
    int used_rects;
private:
    int holes_;
    ll area_;
};

class Solver
{
public:
    Solver(vector<int> a, vector<int> b)
        : n(a.size()), rects(n)
    {
        rep(i, n)
            rects[i] = Rect(a[i], b[i], i);

        ween_creater = WeenCreater(rects);
    }

    int ween_rect;
    int exp_holes;
    int used_rects;
    void solve(int weens)
    {
        ll waste_len = 0;
        ween_rect = 0;
        exp_holes = 1;
        used_rects = 0;
        map<int, int> holes_log;

        vector<bool> used(n);
        Pos wp(-5000, 0);
        rep(ween_i, weens)
        {
            Ween ween;
            if (!ween_creater.create(ween, used))
                break;

//             dump(ween.rects().size());
//             for (auto& r : ween.rects())
//                 fprintf(stderr, "(%4d, %4d)\n", r.long_len(), r.short_len());
//             cerr << endl;

            ween_rect += ween.rects().size();
            waste_len += ween.rects()[0].long_len() + ween.rects()[2].short_len();
            exp_holes += ween.rects().size() - 3;
            holes_log[ween.rects().size() - 3]++;
            used_rects += ween.rects().size();

//             dump(ween.rects().size());

            ween.set_pos(wp);
            for (auto& rect : ween.rects())
            {
                assert(used[rect.index()]);
                assert(check_pos(rect.pos()));
                rects[rect.index()] = rect;
            }

            assert(wp.y <= MAX_RANGE);
            wp += Pos(-5000, 0);
            if (wp.x < -5000 * sqrt(weens))
                wp = Pos(-5000, wp.y + 5000);
        }
#if 0
        int num_ween = 0;
        int sum_holes = 0;
        for (auto& it : holes_log)
        {
            fprintf(stderr, "%3d: %2d\n", it.first, it.second);
            sum_holes += it.first * it.second;
            num_ween += it.second;
        }
        dump(num_ween);
        dump(sum_holes);
        dump(double(sum_holes) / num_ween);
#endif

        vector<int> len_order;
        rep(i, n)
            if (!used[i])
                len_order.push_back(i);
        sort(all(len_order), [&](int i, int j) { return max(rects[i].width(), rects[i].height()) < max(rects[j].width(), rects[j].height()); });
        reverse(all(len_order));
        make_big_square(len_order);

#if 1
#if 0
        for (int i : len_order)
        {
            auto& r = rects[i];
            r.set_pos(Pos(0, 0));
        }
#endif

        dwidth += waste_len / 4;
        dheight += waste_len / 4;
        darea = dwidth * dheight;
#endif
    }

    ll dwidth, dheight;
    ll darea;
    void make_big_square(vector<int> use_rects)
    {
        vector<int> sides[4];
        int sum_len[4] = {};
        for (int rect_i : use_rects)
        {
            int len = max(rects[rect_i].width(), rects[rect_i].height());
            int k = min_element(sum_len, sum_len + 4) - sum_len;
            sum_len[k] += len;
            sides[k].push_back(rect_i);
        }

        pint order[4];
        rep(i, 4)
            order[i] = pint(sum_len[i], i);
        sort(order, order + 4, greater<pint>());

        const int width = order[1].first;
        const int height = order[3].first;
        dwidth = width, dheight = height;

        // 0 bottom
        {
            int x = 0;
            for (int ri : sides[order[0].second])
            {
                int high = max(rects[ri].width(), rects[ri].height());
                int low = min(rects[ri].width(), rects[ri].height());
                if (rects[ri].width() < rects[ri].height())
                    rects[ri].rotate();
                rects[ri].set_pos(Pos(x, -low));

                x += high;
            }
        }

        // 1 top
        {
            int x = 0;
            for (int ri : sides[order[1].second])
            {
                int high = max(rects[ri].width(), rects[ri].height());
                int low = min(rects[ri].width(), rects[ri].height());
                if (rects[ri].width() < rects[ri].height())
                    rects[ri].rotate();
                rects[ri].set_pos(Pos(x, height));

                x += high;
            }
        }

        // 2 left
        {
            int y = 0;
            for (int ri : sides[order[2].second])
            {
                int high = max(rects[ri].width(), rects[ri].height());
                int low = min(rects[ri].width(), rects[ri].height());
                if (rects[ri].width() > rects[ri].height())
                    rects[ri].rotate();
                rects[ri].set_pos(Pos(-low, y));

                y += high;
            }
        }

        // 3 right
        {
            int y = 0;
            for (int ri : sides[order[3].second])
            {
                int high = max(rects[ri].width(), rects[ri].height());
                int low = min(rects[ri].width(), rects[ri].height());
                if (rects[ri].width() > rects[ri].height())
                    rects[ri].rotate();
                rects[ri].set_pos(Pos(width, y));

                y += high;
            }
        }
    }

    Score eval()
    {
#ifndef NDEBUG
        for (auto& rect : rects)
            assert(max(abs(rect.pos().x), abs(rect.pos().y)) <= MAX_RANGE);
#endif

        vector<int> ux, uy;
        ux.push_back(LOWER - 10);
        ux.push_back(UPPER + 10);
        uy.push_back(LOWER - 10);
        uy.push_back(UPPER + 10);
        for (auto& rect : rects)
        {
            ux.push_back(rect.pos().x);
            ux.push_back(rect.pos().x + rect.width());
            uy.push_back(rect.pos().y);
            uy.push_back(rect.pos().y + rect.height());
        }
        uniq(ux);
        uniq(uy);

        const int UNVISIT = -1;
        const int RECT = -1919;
        static int c[2048][2048];
        clr(c, -1);
        for (auto& rect : rects)
        {
            int sx = lower_bound(all(ux), rect.pos().x) - ux.begin();
            int ex = lower_bound(all(ux), rect.pos().x + rect.width()) - ux.begin();
            int sy = lower_bound(all(uy), rect.pos().y) - uy.begin();
            int ey = lower_bound(all(uy), rect.pos().y + rect.height()) - uy.begin();

            for (int y = sy; y < ey; ++y)
                for (int x = sx; x < ex; ++x)
                    c[y][x] = RECT;
        }

        int holes = 0;
        queue<pint> q;
        rep(y, uy.size()) rep(x, ux.size())
        {
            if (c[y][x] == UNVISIT)
            {
                c[y][x] = holes;
                q.push(pint(x, y));
                while (!q.empty())
                {
                    int cx = q.front().first, cy = q.front().second;
                    q.pop();

                    rep(dir, 4)
                    {
                        int nx = cx + DX[dir], ny = cy + DY[dir];
                        if (in_rect(nx, ny, ux.size(), uy.size()) && c[ny][nx] == UNVISIT)
                        {
                            c[ny][nx] = holes;
                            q.push(pint(nx, ny));
                        }
                    }
                }
                ++holes;
            }
        }
        --holes;

        ll area = 0;
        rep(y, uy.size()) rep(x, ux.size())
        {
            if (c[y][x] > 0)
                area += (ll)(uy[y + 1] - uy[y]) * (ux[x + 1] - ux[x]);
        }
//         return Score(holes, area);
        auto s = Score(holes, area);
        s.darea = darea;
        s.dscore = darea * s.holes() * s.holes();
        s.hole_cost =  double(ween_rect) / (s.holes() - 1);
        s.exp_holes = exp_holes;
        s.used_rects = used_rects;
        return s;
    }

    vector<int> make_result()
    {
        vector<int> res;
        rep(i, n)
        {
            auto& r = rects[i];
            res.push_back(r.pos().x);
            res.push_back(r.pos().y);
            res.push_back(r.rotated());
        }
        return res;
    }

private:
    const int n;
    vector<Rect> rects;

    WeenCreater ween_creater;
};


void test(vector<int> a, vector<int> b)
{
    const int n = a.size();
    vector<int> len;
    rep(i, n)
//         len.push_back(max(a[i], b[i]));
        len.push_back(sqrt(a[i] * a[i] + b[i] * b[i]));
    sort(all(len));

    Score best(0, 0);
    rep(i, n)
    {
        ll peri = 0;
        int holes = 1 + i / 1.2;
//         for (int j = 0; j < i; ++j)
//             peri += len[j] / 2;
        for (int j = i; j < n; ++j)
            peri += len[j];

        ll side = peri / 4;
        ll area = side * side;
        Score s(holes, area);
        if (s.score() > best.score())
        {
            best = s;
            fprintf(stderr, "%3d %10lld %13lld\n", holes, area, s.score());
        }
    }
}
void test2(vector<int> a, vector<int> b)
{
    const int n = a.size();
    vector<int> len;
    rep(i, n)
        len.push_back(max(a[i], b[i]));
//         len.push_back(sqrt(a[i] * a[i] + b[i] * b[i]));
    sort(all(len));
    if (n & 1)
        len.push_back(0);

    Score best(0, 0);
    ll peri = accumulate(all(len), 0LL);
    int holes = 1;
    for (int i = 0; i < n; i += 2)
    {
        ll side = peri / 4;
        ll area = side * side;
        Score s(holes, area);
        if (s.score() > best.score())
        {
            best = s;
            fprintf(stderr, "%3d %10lld %13lld\n", holes, area, s.score());
        }

        peri -= len[i + 1];
        ++holes;
    }
}

int arg_weens = 10;
class RectanglesAndHoles
{
public:
    vector<int> place(vector<int> a, vector<int> b)
    {
        const double G_TLE = 9.8 * 1000;
        Timer g_timer;
        g_timer.start();

        Score prev(0, 0);
        const int n = a.size();
        Score best_score(0, 0);
        vector<int> best_res;
        rep(weens, max(1, n / 4 + 1))
//         rep(weens, n + 1)
//         int weens = arg_weens;
        {
//             if (g_timer.get_elapsed() > G_TLE)
//                 break;

            Solver solver(a, b);
            solver.solve(weens);
            Score score = solver.eval();
            if (score.score() > best_score.score())
            {
                best_score = score;
                best_res = solver.make_result();
//                                 fprintf(stderr, "%3d: %3d, %13lld %16lld !\n", weens, score.holes(), score.area(), score.score());
//                                 dump(score.hole_cost);
                fprintf(stderr, "%3d(%3d): %3d[%3d] %+3d %.3f, %11lld[%11lld] %16lld[%16lld] !\n", weens, score.used_rects, score.holes(), score.exp_holes, score.holes() - prev.holes(), score.hole_cost, score.area(),score.darea,  score.score(), score.dscore);
            }
            else
            {
                if (score.score() != prev.score())
//                                         fprintf(stderr, "%3d: %3d, %11lld %16lld\n", weens, score.holes(), score.area(), score.score());
                fprintf(stderr, "%3d(%3d): %3d[%3d] %+3d %.3f, %11lld[%11lld] %16lld[%16lld]\n", weens, score.used_rects, score.holes(), score.exp_holes, score.holes() - prev.holes(), score.hole_cost, score.area(),score.darea,  score.score(), score.dscore);
            }
                prev = score;
        }

        test2(a, b);

        return best_res;
    }
};


#ifdef LOCAL
int main(int argc, char** argv)
{
    if (argc > 1)
        arg_weens = atoi(argv[1]);

    int n;
    cin >> n;
    vector<int> a(n), b(n);
    input(a, n);
    input(b, n);

    vector<int> res = RectanglesAndHoles().place(a, b);
    for (int r : res)
        cout << r << endl;
}
#endif

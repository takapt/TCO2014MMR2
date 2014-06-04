#define NDEBUG

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

template <typename T> void print(T a, int n, const string& split = " ") { for (int i = 0; i < n; i++) { cerr << a[i]; if (i + 1 != n) cerr << split; } cerr << endl; }
template <typename T> void print2d(T a, int w, int h, int width = -1, int br = 0) { for (int i = 0; i < h; ++i) { for (int j = 0; j < w; ++j) { if (width != -1) cerr.width(width); cerr << a[i][j] << ' '; } cerr << endl; } while (br--) cerr << endl; }
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

    bool intersect(const Rect& other) const
    {
        return pos().x < other.pos().x + other.width() && other.pos().x < pos().x + width() &&
        pos().y < other.pos().y + other.height() && other.pos().y < pos().y + height();
    }

private:
    int w, h;
    int index_;
    bool rotated_;
    Pos pos_;
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
private:
    int holes_;
    ll area_;
};

class Solver
{
public:
    Solver(vector<int> a, vector<int> b)
        : n(a.size()), rects(n), used(n)
    {
        rep(i, n)
            rects[i] = Rect(a[i], b[i], i);
    }

    void solve(int num_teeth)
    {
        make_teeth(num_teeth);

        vector<int> sq;
        rep(i, n)
            if (!used[i])
                sq.push_back(i);
        make_big_square(sq);
    }

    void make_teeth(int num_teeth)
    {
        if (num_teeth == 0)
            return;

        vector<int> long_order;
        rep(i, n)
            long_order.push_back(i);
        sort(all(long_order), [&](int i, int j) { return rects[i].long_len() < rects[j].long_len(); });

        vector<int> use_rects;
        rep(i, num_teeth)
        {
            assert(!used[long_order[i]]);
            use_rects.push_back(long_order[i]);
        }
        sort(all(use_rects), [&](int i, int j) { return rects[i].long_len() < rects[j].long_len(); });

        int need_height = 0;
        int fix_width = -1;
        for (int i : use_rects)
        {
            auto& r = rects[i];
            r.width_to_long();
            upmax(fix_width, r.width() + 1);
            need_height += r.height();
        }
        need_height -= rects[use_rects.back()].height();

        const int start_y = 0;
        for (int i = 0, j = sz(use_rects) - 1, y = start_y; i <= j; ++i, --j)
        {
            {
                auto& a = rects[use_rects[i]];
                assert(!used[a.index()]);
                used[a.index()] = true;
                a.set_pos(Pos(fix_width - a.width(), y));
                y += a.height();
            }

            if (i < j)
            {
                auto& b = rects[use_rects[j]];
                assert(!used[b.index()]);
                used[b.index()] = true;
                b.set_pos(Pos(0, y));
                y += b.height();
            }
        }

        for (int i = num_teeth, y = 0; i < n && y < need_height; ++i)
        {
            auto& r = rects[long_order[i]];
            assert(!used[r.index()]);
            used[r.index()] = true;

            r.height_to_long();
            r.set_pos(Pos(fix_width, y));
            y += r.height();
        }
    }

    void make_big_square(vector<int> use_rects)
    {
        sort(all(use_rects), [&](int i, int j) { return rects[i].long_len() < rects[j].long_len(); });
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

        rep(j, n) rep(i, j)
            if (rects[i].intersect(rects[j]))
                return Score(0, 0);

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
    vector<bool> used;
};


class RectanglesAndHoles
{
public:
    vector<int> place(vector<int> a, vector<int> b)
    {
        const double G_TLE = 9.7 * 1000;
        Timer g_timer;
        g_timer.start();

        Score prev(0, 0);
        const int n = a.size();
        Score best_score(0, 0);
        vector<int> best_res;
//         rep(teeth, n)
        for (int teeth = n / 3; teeth < n; ++teeth)
        {
//             if (g_timer.get_elapsed() > G_TLE)
//                 break;

            Solver solver(a, b);
            solver.solve(teeth);
            Score score = solver.eval();
            if (teeth > n / 10 && score.score() == 0)
            {
                dump(teeth);
                break;
            }

            if (score.score() > best_score.score())
            {
                best_score = score;
                best_res = solver.make_result();
                                fprintf(stderr, "%3d: %3d, %13lld %16lld !\n", teeth, score.holes(), score.area(), score.score());
            }
            else
            {
                if (score.score() != prev.score())
                                        fprintf(stderr, "%3d: %3d, %13lld %16lld\n", teeth, score.holes(), score.area(), score.score());
            }
                prev = score;
        }

        return best_res;
    }
};


#ifdef LOCAL
int main(int argc, char** argv)
{
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

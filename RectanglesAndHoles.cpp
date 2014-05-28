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


class InputRect
{
public:
    InputRect(int a, int b, int index)
        : w(a), h(b), index_(index), rotated_(false)
    {
    }
    InputRect()
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

private:
    int w, h;
    int index_;
    bool rotated_;
};



class RectanglesAndHoles
{
public:
    vector<int> place(vector<int> a, vector<int> b)
    {
        const double G_TLE = 10 * 1000;
        Timer g_timer;
        g_timer.start();

        this->a = a;
        this->b = b;
        n = a.size();

        map<int, int> freq;
        rep(i, n)
            ++freq[max(a[i], b[i]) / 100];
        for (auto it : freq)
            fprintf(stderr, "%3d: %d\n", it.first, it.second);

        vector<int> len_order = sorted_by_len();

        vector<int> best_res;
        ll best_score = -114514;
        pair<ll, ll> best_ha;
        for (int tekito = 0; tekito <= n; tekito += 4)
        {
            if (g_timer.get_elapsed() > 9.8 * 1000)
                break;

            clr(rotate, -1);

            vector<int> square, hole;
            rep(i, n)
            {
                if (i < tekito)
                    hole.push_back(len_order[i]);
                else
                    square.push_back(len_order[i]);
            }

            reverse(all(square));
            make_big_square(square);
            make_tekito_holes(hole);

            assert(count(rotate, rotate + n, -1) == 0);
            vector<int> res;
            rep(i, n)
            {
                res.push_back(bottomleft[i].first);
                res.push_back(bottomleft[i].second);
                res.push_back(rotate[i]);
            }

            pair<ll, ll> ha = eval();
            ll score = ha.first * ha.first * ha.second;
            if (score > best_score)
            {
                best_score = score;
                best_ha = ha;
                best_res = make_result();

                fprintf(stderr, "%4d: %3lld, %11lld -> %lld\n", tekito, ha.first, ha.second, score);
            }
        }
        return best_res;
    }

    vector<int> make_result()
    {
#ifndef NDEBUG
        assert(count(rotate, rotate + n, -1) == 0);
        rep(i, n)
            assert(max(abs(bottomleft[i].first), abs(bottomleft[i].second)) <= MAX_RANGE);
#endif

        vector<int> res;
        rep(i, n)
        {
            res.push_back(bottomleft[i].first);
            res.push_back(bottomleft[i].second);
            res.push_back(rotate[i]);
        }
        return res;
    }

    void make_big_square(vector<int> use_rects)
    {
        vector<int> sides[4];
        int sum_len[4] = {};
        for (int rect_i : use_rects)
        {
            int len = max(a[rect_i], b[rect_i]);
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
                int high = max(a[ri], b[ri]);
                int low = min(a[ri], b[ri]);
                assert(rotate[ri] == -1);
                rotate[ri] = a[ri] > b[ri] ? 0 : 1;
                bottomleft[ri] = pint(x, -low);

                x += high;
            }
        }

        // 1 top
        {
            int x = 0;
            for (int ri : sides[order[1].second])
            {
                int high = max(a[ri], b[ri]);
                int low = min(a[ri], b[ri]);
                assert(rotate[ri] == -1);
                rotate[ri] = a[ri] > b[ri] ? 0 : 1;
                bottomleft[ri] = pint(x, height);

                x += high;
            }
        }

        // 2 left
        {
            int y = 0;
            for (int ri : sides[order[2].second])
            {
                int high = max(a[ri], b[ri]);
                int low = min(a[ri], b[ri]);
                assert(rotate[ri] == -1);
                rotate[ri] = a[ri] > b[ri] ? 1 : 0;
                bottomleft[ri] = pint(-low, y);

                y += high;
            }
        }

        // 3 right
        {
            int y = 0;
            for (int ri : sides[order[3].second])
            {
                int high = max(a[ri], b[ri]);
                int low = min(a[ri], b[ri]);
                assert(rotate[ri] == -1);
                rotate[ri] = a[ri] > b[ri] ? 1 : 0;
                bottomleft[ri] = pint(width, y);

                y += high;
            }
        }
    }

    void make_tekito_holes(vector<int> use_rects)
    {
        const int base_x = 200 * 1000;
        int bx = base_x, by = 0;
        for (int ui = 0; ui + 4 <= use_rects.size(); ui += 4)
        {
            rep(i, 4)
                rotate[use_rects[ui + i]] = 0;
            int bottom = use_rects[ui], left = use_rects[ui + 1], top = use_rects[ui + 2], right = use_rects[ui + 3];
            bottomleft[bottom] = pint(bx, by - b[bottom]);
            bottomleft[left] = pint(bx - a[left], by);
            bottomleft[top] = pint(bx, by + min(b[left], b[right]));
            if (a[bottom] > a[top])
            {
                bottomleft[right] = pint(bx + a[top], by);
            }
            else
            {
                bottomleft[right] = pint(bx + a[bottom], by + min(b[left], b[right]) - b[right]);
            }

            const int bd = 5 * 1000;
            bx += bd;
            if (bx >= base_x + 30 * bd)
            {
                bx = base_x;
                by += bd;
            }
        }

        rep(i, use_rects.size() % 4)
        {
            int k = use_rects[use_rects.size() - 1 - i];
            rotate[k] = 0;
            bottomleft[k] = pint(-2000, 2000 * i);
        }
    }

    vector<int> sorted_by_len()
    {
        vector<pint> v;
        rep(i, n)
            v.push_back(pint(max(a[i], b[i]), i));
        sort(all(v));

        vector<int> res;
        rep(i, n)
            res.push_back(v[i].second);
        return res;
    }

    pair<ll, ll> eval()
    {
#ifndef NDEBUG
        assert(count(rotate, rotate + n, -1) == 0);
        rep(i, n)
            assert(max(abs(bottomleft[i].first), abs(bottomleft[i].second)) <= MAX_RANGE);
#endif

        vector<int> ux, uy;
        ux.push_back(LOWER - 10);
        ux.push_back(UPPER + 10);
        uy.push_back(LOWER - 10);
        uy.push_back(UPPER + 10);
        rep(i, n)
        {
            int w = rotate[i] == 0 ? a[i] : b[i];
            int h = rotate[i] == 0 ? b[i] : a[i];
            ux.push_back(bottomleft[i].first);
            ux.push_back(bottomleft[i].first + w);
            uy.push_back(bottomleft[i].second);
            uy.push_back(bottomleft[i].second + h);
        }
        uniq(ux);
        uniq(uy);

        const int UNVISIT = -1;
        const int RECT = -1919;
        static int c[2048][2048];
        clr(c, -1);
        rep(i, n)
        {
            int w = rotate[i] == 0 ? a[i] : b[i];
            int h = rotate[i] == 0 ? b[i] : a[i];
            int sx = lower_bound(all(ux), bottomleft[i].first) - ux.begin();
            int ex = lower_bound(all(ux), bottomleft[i].first + w) - ux.begin();
            int sy = lower_bound(all(uy), bottomleft[i].second) - uy.begin();
            int ey = lower_bound(all(uy), bottomleft[i].second + h) - uy.begin();

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
        return pair<ll, ll>(holes, area);
    }

    int n;
    vector<int> a, b;

    pint bottomleft[1024];
    int rotate[1024];
};


#ifdef LOCAL
int main()
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

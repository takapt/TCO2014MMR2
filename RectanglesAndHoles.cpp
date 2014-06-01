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

        rep(i, 4)
            outer[i].set_pos(outer[i].pos() + move);
        for (auto& rect : inner)
            rect.set_pos(rect.pos() + move);
    }

    vector<Rect> rects() const
    {
        vector<Rect> res;
        rep(i, 4)
            res.push_back(outer[i]);
        for (auto& rect : inner)
            res.push_back(rect);

        for (auto& rect : res)
            rect.set_pos(pos_ + rect.pos());
        return res;
    }

private:
    // positions are relative to pos_
    Rect outer[4];
    vector<Rect> inner;

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
        outer_cand.resize(rects.size());
        rep(i, outer_cand.size())
            outer_cand[i] = i;
        sort(all(outer_cand), [&](int i, int j) { return rects[i].long_len() < rects[j].long_len(); });

        // TODO: (long, index), (short, index)の両方をinner_candに突っ込んで決めたほうがいいか？
        rep(i, rects.size())
            inner_cand.insert(pint(rects[i].long_len(), i));
    }

    bool create(Ween& ween, vector<bool>& used)
    {
        ween = Ween();

        //
        rep(i, 4)
        {
            while (!outer_cand.empty() && used[rects[outer_cand.back()].index()])
                outer_cand.pop_back();
            if (outer_cand.empty())
                return false;

            ween.outer[i] = rects[outer_cand.back()];
            outer_cand.pop_back();
        }
        rep(i, 4)
            used[ween.outer[i].index()] = true;

        // 
        Rect& bottom = ween.outer[0];
        bottom.width_to_long();

        Rect& top = ween.outer[1];
        top.width_to_long();

        Rect& left = ween.outer[2];
        left.height_to_long();

        Rect& right = ween.outer[3];
        right.height_to_long();

        bottom.set_pos(Pos(0, -bottom.height()));
        top.set_pos(Pos(0, right.height()));
        left.set_pos(Pos(-left.width(), 0));
        create_recur(right, top.width() + right.width() + 1, right.height(), MoveDir::X, ween, used);

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
                return;

            next_ween_rect.width_to_long();
            ween_rect.set_pos(Pos(next_ween_rect.width(), 0));
            used[ween_rect.index()] = true;
            ween.inner.push_back(ween_rect);

            create_recur(next_ween_rect, next_ween_rect.width(), height, MoveDir::Y, ween, used);
        }
        else if (move_dir == Y)
        {
            const int next_height_range = height - ween_rect.height() - 1;
            assert(next_height_range > 0);
            Rect next_ween_rect;
            if (!select_ween_rect(next_height_range, width - 2, next_ween_rect, used))
                return;

            next_ween_rect.height_to_long();
            ween_rect.set_pos(Pos(0, next_ween_rect.height()));
            used[ween_rect.index()] = true;
            ween.inner.push_back(ween_rect);

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

    void solve()
    {
        vector<bool> used(n);
        Ween ween;
        dump(ween_creater.create(ween, used));
        ween.set_pos(Pos(-5000, -5000));
        for (auto& rect : ween.rects())
            rects[rect.index()] = rect;

        vector<int> len_order;
        rep(i, n)
            if (!used[i])
                len_order.push_back(i);
        sort_by_long_side(len_order);
        reverse(all(len_order));
        make_big_square(len_order);
    }

    bool check_rect_indices(const vector<int>& rect_i)
    {
        if (rect_i.empty())
            return true;
        auto v = rect_i;
        uniq(v);
        return v.size() == rect_i.size() && 0 <= v.front() && v.back() < n;
    }

    void sort_by_long_side(vector<int>& order)
    {
        assert(check_rect_indices(order));
        sort(all(order), [&](int i, int j) { return max(rects[i].width(), rects[i].height()) < max(rects[j].width(), rects[j].height()); });
    }

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


class RectanglesAndHoles
{
public:
    vector<int> place(vector<int> a, vector<int> b)
    {
        Solver solver(a, b);
        solver.solve();
        return solver.make_result();
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

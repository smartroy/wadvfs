#include <queue>
#include <iostream>
using namespace std;

int main()
{
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> q;

    q.push({6, 7});
    q.push({6, 6});
    q.push({6, 2});
    q.push({6, 3});
    q.push({6, 4});
    while (!q.empty())
    {
        auto p = q.top();
        q.pop();
        std::cout << p.first << " " << p.second << std::endl;
    }
}
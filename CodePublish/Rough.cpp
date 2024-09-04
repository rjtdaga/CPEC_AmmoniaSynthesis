#include <iostream>
using namespace std;

enum seasons { spring = 34, summer = 34.1, autumn = 9, winter = 32};

int main() {

    seasons s;

    s = summer;
    if (summer == spring)
      cout << "Yes" << endl;
    else
       cout << "No" << endl;

    return 0;
}
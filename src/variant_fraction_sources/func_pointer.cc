#include <cstdio>

int add(int a,int b)
{

  return a+b;

}

int sub(int a,int b)
{

  return a-b;

}

//関数ポインタを使っている
int calc(int (*func)(int,int),int a,int b)
{

  return func(a,b);

}

int main(void)
{

  int x = 10;

  int y = -4;

  int ans1 = calc(add,x,y);
  //足し算をしている
  int ans2 = calc(sub,x,y);
  //引き算をしている
  printf("%d,%d\n",ans1,ans2);
  //6,14が出力される
  return 0;

}

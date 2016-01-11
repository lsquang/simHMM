#include "COneSite.h"

bool operator == (const onesite &x, const onesite &y)
{
  return (x.pos == y.pos && x.al1 == y.al1 && x.al2 == y.al2);
}

bool operator < (const onesite &x, const onesite &y)
{
  if (x.pos != y.pos) return (x.pos < y.pos);
  if (x.al1 != y.al1) return (x.al1 < y.al1);
  return (x.al2 < y.al2);
}


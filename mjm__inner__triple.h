

// no guard etc deisned for inner class includsion

class _triple
{
typedef D Tn;
public:
// don't bother init... doh 
_triple() {}
_triple(const Tn & x, const Tn &y, const Tn & z): m_x(x),m_y(y),m_z(z) {}
_triple(const Tn & x, const Tn &y): m_x(x),m_y(y),m_z(0) {}
const Tn & x() const { return m_x;}
const Tn & y() const { return m_y;}
const Tn & z() const { return m_z;}
//Tn & operator[](const IdxTy i) { return *(((Tn*)this)+i); } 
Tn & operator[](const IdxTy i) { return (&m_x)[i]; } 
const Tn & operator[](const IdxTy i)const { return (&m_x)[i]; } 
_triple & operator=(const D &  x) {m_x=x; m_y=x; m_z=x;  return *this ; } 
bool operator!=(const D &  x)const 
{
if (m_x!=x) return true;
if (m_y!=x) return true;
if (m_z!=x) return true;
 return false ; 
} 

_triple& move(const Tn & x, const Tn &y, const Tn & z)
{
m_x+=x; m_y+=y; m_z+=z; 
return *this;
} // move

_triple& scale(const D & s)// const
{
_triple& x=(*this);
x.m_x*=s; x.m_y*=s; x.m_z*=s;
return x; 
} 
_triple operator*(const D & s) const
{
_triple x=(*this);
x.m_x*=s; x.m_y*=s; x.m_z*=s;
return x; 
} 
_triple operator/(const D & s) const
{
_triple x=(*this);
x.m_x/=s; x.m_y/=s; x.m_z/=s;
return x; 
} 


_triple operator+(const _triple & s) const
{
_triple x=(*this);
x.m_x+=s.m_x; x.m_y+=s.m_y; x.m_z+=s.m_z;
return x; 
} 
_triple operator/(const _triple & s) const
{
_triple x=(*this);
x.m_x/=s.m_x; x.m_y/=s.m_y; x.m_z/=s.m_z;
return x; 
} 


_triple operator-(const _triple & s) const
{
_triple x=(*this);
x.m_x-=s.m_x; x.m_y-=s.m_y; x.m_z-=s.m_z;
return x; 
} 
_triple operator*(const _triple & s) const
{
_triple x=(*this);
x.m_x*=s.m_x; x.m_y*=s.m_y; x.m_z*=s.m_z;
return x; 
} 


_triple& operator+=(const _triple & s) 
{
_triple& x=(*this);
x.m_x+=s.m_x; x.m_y+=s.m_y; x.m_z+=s.m_z;
return x; 
} 
template <class Ta> _triple& operator+=(const Ta & s) 
{
_triple& x=(*this);
x.m_x+=s[0]; x.m_y+=s[1]; x.m_z+=s[2];
return x; 
} 



_triple& operator-=(const _triple & s) 
{
_triple& x=(*this);
x.m_x-=s.m_x; x.m_y-=s.m_y; x.m_z-=s.m_z;
return x; 
} 



Tn dot(const _triple & s) const
{
Tn x=x.m_x*s.m_x+ x.m_y*s.m_y+x.m_z*s.m_z;
return x; 
} 
_triple cross(const _triple & s) const
{
_triple x; // (*this);
x.m_x=m_y*s.m_z-s.m_y*m_z;// 
x.m_y=-m_x*s.m_z+s.m_x*m_z;// 
x.m_z=m_x*s.m_y-s.m_x*m_y;// 
return x; 
} 

_triple& unit() 
{
_triple& x=(*this);
const Tn d=sqrt(r2());
if (d==0) return x; 
// should make inverse and only div once... 
//x.m_x/=d; x.m_y/=d; x.m_z/d;
x=x*(1.0/d);
return x; 
} 

Tn r2() const { return m_x*m_x+m_y*m_y+m_z*m_z; } 
StrTy dump() const
{ Ss ss; dump(ss); return ss.str(); } 
void dump (OsTy& os) const
{ os<<MMPR4(m_x,m_y,m_z,r2()); } 

Tn  m_x,m_y,m_z;
}; // _triple

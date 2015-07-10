#ifndef Box_h
#define Box_h

#include "containers/vec.h"

using namespace LMT;
using namespace std;

template<unsigned dim_,class TT_> class Box
{

    public:

        Box()
        {
            isCornersSet=false;
            if ( dim_==2 )
            {
                base[0]=Vec<TT_,2> ( 1.,0. );
                base[1]=Vec<TT_,2> ( 0.,1. );
            }
            else if ( dim_==3 )
            {
                base[0]=Vec<TT_,3> ( 1.,0.,0. );
                base[1]=Vec<TT_,3> ( 0.,1.,0. );
                base[2]=Vec<TT_,3> ( 0.,0.,1. );
            }
        }

        void setCorners ( Vec<Vec<TT_,dim_>,2>& corners_ )
        {
            corners=corners_;
            isCornersSet=true;
        }

        const Vec<Vec<TT_,dim_>,2>& getCorners() const
        {
            if ( isCornersSet )
                return corners;
            else assert ( 0 /*Les coins de la boites ne sont pas définis */ );
        }

        /** \ingroup  Maillages
        \brief Recherche si un pt est a l'interieur d'une boite d�finie par ses deux extr�mit�s et par une base dont le dernier vecteur est la normale. Une proc�dure est �crite en 2d et une autre en 3d

        Pour le cas 2d, on cherche � d�terminer si le pt M est dans le segment [AB]. Pour cela, on calcul les distances MA et MB puis on v�rifie si \f$ \frac{\vec{MA}.\vec{MB}}{MA MB} =-1 \f$ � epsilon pr�s.

        Pour le cas 3d, connaissant la base dont le dernier vecteur est la normale � la surface et les deux autres sont parall�les aux cot�s, on cherche � savoir si le pt M est dans le rectangle d�fini par les deux points extr�mit�s de la boite. On calcule donc :
        \f$ (\vec{MA}.\vec{n})x (\vec{MB}.\vec{n}) <=0 \f$
        \f$ (\vec{MA}.\vec{t1})x (\vec{MB}.\vec{t1}) <=0 \f$
        \f$ (\vec{MA}.\vec{t2})x (\vec{MB}.\vec{t2}) <=0 \f$
        Attention il faut n�cessairement que t1 et t2 soient parall�les aux cot�s du rectangle (� construire lorsque l'on cr�e la base).
        */
//  bool pt_in_box(const Vec<TT_,2> &pt, Vec<Vec<TT_,2>,2> &base, const double eps=1e-6) const {
//      TT_ la = length(pt-corners[0]);
//      TT_ lb = length(pt-corners[1]);
//      bool flag=false;
//      if(la<=eps or lb<=eps or abs(dot(pt-corners[0],pt-corners[1])+la*lb)<=eps)
//              flag=true;
//      return flag;
//  }
//  bool pt_in_box(const Vec<TT_,3> &pt, Vec<Vec<TT_,3>,3> &base, const double eps=1e-6) const {
//      bool flag=false;
//      if( (dot(pt-corners[0],base[2])*dot(pt-corners[1],base[2])<=eps) and (dot(pt-corners[0],base[0])*dot(pt-corners[1],base[0])<=eps) and (dot(pt-corners[0],base[1])*dot(pt-corners[1],base[1])<=eps) )
//          flag=true;
//      return flag;
//  }

        /** \ingroup  Maillages
        \brief Recherche si un pt est a l'interieur d'une boite d�finie par ses deux extr�mit�s (rectangle en 2d et parall�l�pip�de en 3d)

        Pour le cas 2d ou 3d, on cherche � savoir si le pt M est dans le rectangle ou (parall�l�pip�de) d�fini par les deux points extr�mit�s de la boite A et B. On calcule donc :
        \f$ (\vec{MA}.\vec{x})x (\vec{MB}.\vec{x}) <=0 \f$
        \f$ (\vec{MA}.\vec{y})x (\vec{MB}.\vec{y}) <=0 \f$
        \f$ (\vec{MA}.\vec{z})x (\vec{MB}.\vec{z}) <=0 \f$
        Attention il faut n�cessairement que t1 et t2 soient parall�les aux cot�s du rectangle (� construire lorsque l'on cr�e la base).
        */
        bool pt_in_box ( const Vec<TT_,2> &pt, const double eps=1e-6 ) const
        {
            bool flag=false;
            Vec<TT_,2> x ( 1.,0. ),y ( 0.,1. );
            if ( ( dot ( pt-corners[0],x ) *dot ( pt-corners[1],x ) <=eps ) and ( dot ( pt-corners[0],y ) *dot ( pt-corners[1],y ) <=eps ) )
                flag=true;
            return flag;
        }
        bool pt_in_box ( const Vec<TT_,3> &pt, const double eps=1e-6 ) const
        {
            bool flag=false;
            Vec<TT_,3> x ( 1.,0.,0. ),y ( 0.,1.,0. ),z ( 0.,0.,1. );
            if ( ( dot ( pt-corners[0],x ) *dot ( pt-corners[1],x ) <=eps ) and ( dot ( pt-corners[0],y ) *dot ( pt-corners[1],y ) <=eps ) and ( dot ( pt-corners[0],z ) *dot ( pt-corners[1],z ) <=eps ) )
                flag=true;
            return flag;
        }

        /** \brief Intersection entre deux boites definies par leurs noeuds extremite. Cette fonction renvoit 1 si les deux boites s'intersectent ou se touchent (� epsilon pr�s), 0 sinon

        Une fonction est �crite pour des boites 2d et une autre pour des boites 3d. Normalement aucune restriction n'est n�cessaire dans l'ordre des points extr�mit�s des boites. Les boites doivent simplement pouvoir �tre repr�sent�es par 2 points et donc align�es sur les axes x, y, z. Cette fonction ne n�cessite en 3d que 3 comparaisons et 2 tests ainsi que l'utilisation de produits scalaires qui sont normalement optimis�s.
        **/
        bool intersection_box ( Box<dim_,TT_> &box2, const double &eps=1e-6 ){
            typedef Vec<TT_,dim_> TV;
            Vec<Vec<TT_,dim_>,2> corners2=box2.getCorners();
            //boite 1
            TV P1 = ( corners[0]+corners[1] ) /2.;
            TV E1 = abs ( corners[1]-corners[0] ) /2.;
            //boite 2
            TV P2 = ( corners2[0]+corners2[1] ) /2.;
            TV E2 = abs ( corners2[1]-corners2[0] ) /2.;

            TV V = P2-P1;
            if (dim_==2)
            {
                TV x ( 1,0 );
                TV y ( 0,1 );
                return ( ( abs ( dot ( V,x ) ) <=dot ( E1,x ) +dot ( E2,x ) +eps ) and ( abs ( dot ( V,y ) ) <=dot ( E1,y ) +dot ( E2,y ) +eps ) ) ;
            }
            else if (dim_==3)
            {
                TV x(1,0,0);
                TV y(0,1,0);
                TV z(0,0,1);
                return ( (abs(dot(V,x))<=dot(E1,x)+dot(E2,x)+eps) and (abs(dot(V,y))<=dot(E1,y)+dot(E2,y)+eps) and (abs(dot(V,z))<=dot(E1,z)+dot(E2,z)+eps) );
            }
        }

            /** \ingroup Maillage_geometrie_inter
            \brief Intersection d'un maillage de peau avec une boite definie par deux points extr�mit� : sort les numéros des éléments
            */
            template<class TM1>
            Vec<unsigned> intersection_mesh_box ( TM1 &m ) const
            {
                //recherche des elements contenus dans la boite
                Vec<unsigned> num_elem;
                for ( unsigned i=0;i<m.elem_list.size();i++ )
                {
                    typename TM1::Pvec G;
                    G = center ( *m.elem_list[i] );
                    if ( pt_in_box ( G ) ==1 )
                        num_elem.push_back ( i );
                }
                return num_elem;
            }

        private:
            Vec<Vec<TT_,dim_>,2> corners;
            bool isCornersSet;
            Vec<Vec<TT_,dim_>,dim_> base;

        };

#endif //Box_h

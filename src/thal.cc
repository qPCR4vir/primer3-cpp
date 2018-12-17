/*
 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2009,2010,
               2011,2012
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
 All rights reserved.

       This file is part of primer3 software suite.

       This software suite is is free software;
       you can redistribute it and/or modify it under the terms
       of the GNU General Public License as published by the Free
       Software Foundation; either version 2 of the License, or (at
       your option) any later version.

       This software is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this software (file gpl-2.0.txt in the source
       distribution); if not, write to the Free Software
       Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <fstream>
#include <strstream>
#include <exception>
#include <map>
#include <vector>

#include "primer3_config/dangle.dh.hpp"
#include "primer3_config/dangle.ds.hpp"
#include "primer3_config/loops.dh.hpp"
#include "primer3_config/loops.ds.hpp"
#include "primer3_config/stack.dh.hpp"
#include "primer3_config/stack.ds.hpp"
#include "primer3_config/stackmm.dh.hpp"
#include "primer3_config/stackmm.dh.hpp"
#include "primer3_config/tetraloop.dh.hpp"
#include "primer3_config/tetraloop.ds.hpp"
#include "primer3_config/triloop.dh.hpp"
#include "primer3_config/triloop.ds.hpp"
#include "primer3_config/tstack.dh.hpp"
#include "primer3_config/tstack_tm_inf.ds.hpp"
#include "primer3_config/tstack2.dh.hpp"
#include "primer3_config/tstack2.ds.hpp"

#if defined(__sun)
#include <ieeefp.h>
#endif
/*#define DEBUG*/   // todo ??

#include "thal.hpp"

extern const double _INFINITY;     // todo ??

//     **************    thal_parameters    ************************


#ifdef INTEGER
# define isFinite(x) (x < _INFINITY / 2)
#else
# define isFinite(x) isfinite(x)
#endif

#define isPositive(x) ((x) > 0 ? (1) : (0))

/*** BEGIN CONSTANTS ***/
# ifdef INTEGER
const double _INFINITY = 999999.0;
# else
#   ifdef INFINITY
const double _INFINITY = INFINITY;
#   else
const double _INFINITY = 1.0 / 0.0;
#   endif
# endif

/* part of calculating salt correction for Tm by SantaLucia et al */
static double
saltCorrectS (double mv, double dv, double dntp)
{
    if(dv<=0) dntp=dv;
    return 0.368*((log((mv+120*(sqrt(fmax(0.0, dv-dntp))))/1000)));
}


/* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
static const int BPI[5][5] =  {
        {0, 0, 0, 1, 0}, /* A, C, G, T, N; */
        {0, 0, 1, 0, 0},
        {0, 1, 0, 0, 0},
        {1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}};

static unsigned char
nt2code(char c)          ///< converts DNA nt char to unsigned char code; 0-A, 1-C, 2-G, 3-T, 4-whatever
{
    switch (c) {
        case 'A': case '0':       return 0;
        case 'C': case '1':       return 1;
        case 'G': case '2':       return 2;
        case 'T': case '3':       return 3;
    }
    return 4;
}

seq nt2code (const std::string& nt)
{
    seq code;
    code.reserve( nt.length() );
    for (char n:nt) code+=nt2code(n);
    return code;
}

static void
reverse(seq& s)                                // not complement !
{
    for (int i = 0, j = s.length()-1; i < j; i++, j--)
    {
        char c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

class thal_parameters::impl
{
    using upis=std::unique_ptr<std::istream>;

    upis dangle_dh,   dangle_ds;
    upis loops_dh,    loops_ds;
    upis stack_dh,    stack_ds;
    upis stackmm_dh,  stackmm_ds;
    upis tetraloop_dh,tetraloop_ds;
    upis triloop_dh,  triloop_ds;
    upis tstack_tm_inf_ds, tstack_dh;
    upis tstack2_dh,  tstack2_ds;

public:
    impl()
    {
        set_defaults( );
    }
    explicit impl(const std::filesystem::path& dirname )
    {
        load(dirname);
    }

    int set_defaults( )
    {
        this->dangle_dh         = std::make_unique<std::istrstream>(dangle_dh_data );
        this->dangle_ds         = std::make_unique<std::istrstream>(dangle_ds_data);
        this->loops_dh          = std::make_unique<std::istrstream>(loops_dh_data );
        this->loops_ds          = std::make_unique<std::istrstream>(loops_ds_data );
        this->stack_dh          = std::make_unique<std::istrstream>(stack_dh_data );
        this->stack_ds          = std::make_unique<std::istrstream>(stack_ds_data );
        this->stackmm_dh        = std::make_unique<std::istrstream>(stackmm_dh_data);
        this->stackmm_ds        = std::make_unique<std::istrstream>(stackmm_ds_data);
        this->tetraloop_dh      = std::make_unique<std::istrstream>(tetraloop_dh_data );
        this->tetraloop_ds      = std::make_unique<std::istrstream>(tetraloop_ds_data );
        this->triloop_dh        = std::make_unique<std::istrstream>(triloop_dh_data );
        this->triloop_ds        = std::make_unique<std::istrstream>(triloop_ds_data );
        this->tstack_tm_inf_ds  = std::make_unique<std::istrstream>(tstack_tm_inf_ds_data );
        this->tstack_dh         = std::make_unique<std::istrstream>(tstack_dh_data );
        this->tstack2_dh        = std::make_unique<std::istrstream>(tstack2_dh_data );
        this->tstack2_ds        = std::make_unique<std::istrstream>(tstack2_ds_data );
    }

    int load(const std::filesystem::path& dirname )
    {
        dangle_dh    = readParamFile(dirname, "dangle.dh");
        dangle_ds    = readParamFile(dirname, "dangle.ds");
        loops_dh     = readParamFile(dirname, "loops.dh");
        loops_ds     = readParamFile(dirname, "loops.ds");
        stack_dh     = readParamFile(dirname, "stack.dh");
        stack_ds     = readParamFile(dirname, "stack.ds");
        stackmm_dh   = readParamFile(dirname, "stackmm.dh");
        stackmm_ds   = readParamFile(dirname, "stackmm.ds");
        tetraloop_dh = readParamFile(dirname, "tetraloop.dh");
        tetraloop_ds = readParamFile(dirname, "tetraloop.ds");
        triloop_dh   = readParamFile(dirname, "triloop.dh");
        triloop_ds   = readParamFile(dirname, "triloop.ds");
        tstack_tm_inf_ds = readParamFile(dirname, "tstack_tm_inf.ds");
        tstack_dh    = readParamFile(dirname, "tstack.dh");
        tstack2_dh   = readParamFile(dirname, "tstack2.dh");
        tstack2_ds   = readParamFile(dirname, "tstack2.ds");

        return 0;
    }

    int  parse_thermodynamic_values( )
    {
        getStack    ();
        getStackint2();
        getDangle   ();
        getLoop     ();
        getTstack   ();
        getTstack2  ();
        getTriloop  ();
        getTetraloop();

        tableStartATS(AT_S, atpS);    /* getting the AT-penalties */
        tableStartATH(AT_H, atpH);

        return 0;
    }

    double  dangleEntropies3   [5][5][5],     ///< thermodynamic paramteres for 3' dangling ends
            dangleEnthalpies3  [5][5][5],     ///< ther params for 3' dangling ends
            dangleEntropies5   [5][5][5],     ///< ther params for 5' dangling ends
            dangleEnthalpies5  [5][5][5],     ///< ther params for 5' dangling ends
            stackEntropies     [5][5][5][5],  ///< ther params for perfect match pairs
            stackEnthalpies    [5][5][5][5],  ///< ther params for perfect match pairs
            stackint2Entropies [5][5][5][5],  ///< ther params for perfect match and internal mm
            stackint2Enthalpies[5][5][5][5],  ///< ther params for perfect match and internal mm
            tstackEntropies    [5][5][5][5],  ///< ther params for terminal mismatches
            tstackEnthalpies   [5][5][5][5],  ///< ther params for terminal mismatches
            tstack2Entropies   [5][5][5][5],  ///< ther params for internal terminal mismatches
            tstack2Enthalpies  [5][5][5][5],  ///< ther params for internal terminal mismatches
            interiorLoopEntropies [30],       ///< interior loop params according to length of the loop
            bulgeLoopEntropies    [30],       ///< bulge loop params according to length of the loop
            hairpinLoopEntropies  [30],       ///< hairpin loop params accordint to length of the loop
            interiorLoopEnthalpies[30],       ///< same as interiorLoopEntropies but values of entropy
            bulgeLoopEnthalpies   [30],       ///< same as bulgeLoopEntropies but values of entropy
            hairpinLoopEnthalpies [30],       ///< same as hairpinLoopEntropies but values of entropy
            atpS            [5][5],           ///< AT penalty
            atpH            [5][5];           ///< AT penalty

    using loop_prmtr = std::map<seq, double> ;   ///< loops parameter as map of uchar sequence code to value
    loop_prmtr  triloopEntropies,      ///< therm penalties for given triloop   seq-s
                triloopEnthalpies,     ///< therm penalties for given triloop   seq-s
                tetraloopEntropies,    ///< therm penalties for given tetraloop seq-s
                tetraloopEnthalpies ;  ///< therm penalties for given tetraloop seq-s

private:
    static upis readParamFile(const std::filesystem::path& dirname, const std::filesystem::path& fname)
    {
        std::ifstream file {(dirname / fname).string()};
        if (!file) {
            throw std::runtime_error( "Trying to read Th parameters file: Unable to open file: "
                                      + (dirname / fname).string());
        }
        // read the file
        auto r = std::make_unique<std::strstream>();
        (*r) << fi.rdbuf();   // this will eats the new-lines ??    //upis res= std::move(r);
        // see https://en.cppreference.com/w/cpp/io/basic_ostream/operator_ltlt
        // http://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
        return r;
    }

    /// These functions are needed as "inf" cannot be read on Windows directly
    static double readDouble(std::istream& istr)
    {
        str::string nmb;
        istr >> nmb ;

        if (nmb == "inf") return _INFINITY;
        return std::stod(nmb);
    }

    /// Reads a line containing 3 doubles, which can be specified as "inf".
    static void readLoop(std::istream& istr , double &v1, double &v2, double &v3)
    {
        int n;         /* skip first number on the line */
        istr >> n ;

        v1 = readDouble(istr);
        v2 = readDouble(istr);
        v3 = readDouble(istr);
    }

    /// Reads a line containing a short string and a double, used for reading a triloop or tetraloop.
    static std::istream& readTLoop(std::istream& istr, loop_prmtr& lprmtr)
    {
        std::string& s; /*tetraloop string has 6 characters*/ /*triloop string has 5 characters*/
        s.reserve(6);
        istr >> s;                            /* read the string */
        lprmtr[ nt2code(s) ] = readDouble(istr);
        return istr;
    }

    void getStack()
    {
        stack_ds->seekg(std::ios_base::beg);
        stack_dh->seekg(std::ios_base::beg);
        for (int i  = 0; i  < 5; ++i )
        for (int ii = 0; ii < 5; ++ii)
        for (int j  = 0; j  < 5; ++j )
        for (int jj = 0; jj < 5; ++jj)
            if (i == 4 || j == 4 || ii == 4 || jj == 4)
            {
                stackEntropies [i][ii][j][jj] = -1.0;
                stackEnthalpies[i][ii][j][jj] = _INFINITY;
            } else
            {
                stackEntropies [i][ii][j][jj] = readDouble(*stack_ds);
                stackEnthalpies[i][ii][j][jj] = readDouble(*stack_dh);
                if (!isFinite(stackEntropies [i][ii][j][jj]) ||
                    !isFinite(stackEnthalpies[i][ii][j][jj]))
                {
                    stackEntropies [i][ii][j][jj] = -1.0;
                    stackEnthalpies[i][ii][j][jj] = _INFINITY;
                }
            }
    }

    void getStackint2()
    {
        stackmm_ds->seekg(std::ios_base::beg);
        stackmm_dh->seekg(std::ios_base::beg);
        for (int i  = 0; i  < 5; ++i )
        for (int ii = 0; ii < 5; ++ii)
        for (int j  = 0; j  < 5; ++j )
        for (int jj = 0; jj < 5; ++jj)
            if (i == 4 || j == 4 || ii == 4 || jj == 4)
            {
                stackint2Entropies [i][ii][j][jj] = -1.0;
                stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
            } else
            {
                stackint2Entropies [i][ii][j][jj] = readDouble(*stackmm_ds );
                stackint2Enthalpies[i][ii][j][jj] = readDouble(*stackmm_dh );
                if (!isFinite(stackint2Entropies [i][ii][j][jj]) ||
                    !isFinite(stackint2Enthalpies[i][ii][j][jj]))
                {
                    stackint2Entropies [i][ii][j][jj] = -1.0;
                    stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
                }
            }
    }

    void getDangle()
    {
        dangle_ds->seekg(std::ios_base::beg);
        dangle_dh->seekg(std::ios_base::beg);
        for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
        for (int k = 0; k < 5; ++k)
        {
            if (i == 4 || j == 4)
            {
                dangleEntropies3 [i][k][j] = -1.0;
                dangleEnthalpies3[i][k][j] = _INFINITY;
            }
            else
            if (k == 4)
              {
                dangleEntropies3 [i][k][j] = -1.0;
                dangleEnthalpies3[i][k][j] = _INFINITY;
              } else
              {
                dangleEntropies3 [i][k][j] = readDouble(*dangle_ds );
                dangleEnthalpies3[i][k][j] = readDouble(*dangle_dh );
                if(!isFinite(dangleEntropies3 [i][k][j]) ||
                   !isFinite(dangleEnthalpies3[i][k][j]))
                {
                    dangleEntropies3 [i][k][j] = -1.0;
                    dangleEnthalpies3[i][k][j] = _INFINITY;
                }
              }
        }

        for (i = 0; i < 5; ++i)
        for (j = 0; j < 5; ++j)
        for (k = 0; k < 5; ++k)
        {
            if (i == 4 || j == 4)
            {
                dangleEntropies5 [i][j][k] = -1.0;
                dangleEnthalpies5[i][j][k] = _INFINITY;
            } else
            if (k == 4)
              {
                dangleEntropies5 [i][j][k] = -1.0;
                dangleEnthalpies5[i][j][k] = _INFINITY;
              } else
              {
                dangleEntropies5 [i][j][k] = readDouble(*tp->dangle_ds );
                dangleEnthalpies5[i][j][k] = readDouble(*tp->dangle_dh);
                if(!isFinite(dangleEntropies5 [i][j][k]) ||
                   !isFinite(dangleEnthalpies5[i][j][k]))
                {
                    dangleEntropies5 [i][j][k] = -1.0;
                    dangleEnthalpies5[i][j][k] = _INFINITY;
                }
              }
        }
    }

    void getLoop()
    {
        loops_ds->seekg(std::ios_base::beg);
        loops_dh->seekg(std::ios_base::beg);
        for (int k = 0; k < 30; ++k)
        {
            readLoop(*loops_ds, &interiorLoopEntropies[k],
                                &bulgeLoopEntropies   [k],
                                &hairpinLoopEntropies [k]);
            readLoop(*loops_ds, &interiorLoopEnthalpies[k],
                                &bulgeLoopEnthalpies   [k],
                                &hairpinLoopEnthalpies [k]);
        }
    }

    void getTstack()
    {
        tstack_tm_inf_ds->seekg(std::ios_base::beg);
        tstack_dh->seekg(std::ios_base::beg);
        for (int i1 = 0; i1 < 5; ++i1)
        for (int i2 = 0; i2 < 5; ++i2)
        for (int j1 = 0; j1 < 5; ++j1)
        for (int j2 = 0; j2 < 5; ++j2)
            if (i1 == 4 || j1 == 4)
            {
                tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
                tstackEntropies [i1][i2][j1][j2] = -1.0;
            } else
              if (i2 == 4 || j2 == 4)
            {
                tstackEntropies [i1][i2][j1][j2] = 0.00000000001;
                tstackEnthalpies[i1][i2][j1][j2] = 0.0;
            } else
            {
                tstackEntropies [i1][i2][j1][j2] = readDouble(*tstack_tm_inf_ds);
                tstackEnthalpies[i1][i2][j1][j2] = readDouble(*tstack_dh);
                if ( !isFinite(tstackEntropies [i1][i2][j1][j2]) ||
                     !isFinite(tstackEnthalpies[i1][i2][j1][j2]) )
                {
                    tstackEntropies [i1][i2][j1][j2] = -1.0;
                    tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
                }
            }
    }

    void getTstack2()
    {
        tstack2_ds->seekg(std::ios_base::beg);
        tstack2_dh->seekg(std::ios_base::beg);
        for (int i1 = 0; i1 < 5; ++i1)
        for (int i2 = 0; i2 < 5; ++i2)
        for (int j1 = 0; j1 < 5; ++j1)
        for (int j2 = 0; j2 < 5; ++j2)
            if (i1 == 4 || j1 == 4)
            {
                tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
                tstack2Entropies [i1][i2][j1][j2] = -1.0;
            } else
            if (i2 == 4 || j2 == 4)
            {
                tstack2Entropies [i1][i2][j1][j2] = 0.00000000001;
                tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
            } else
            {
                tstack2Entropies [i1][i2][j1][j2] = readDouble(*tstack2_ds);
                tstack2Enthalpies[i1][i2][j1][j2] = readDouble(*tstack2_dh);
                if (!isFinite(tstack2Entropies [i1][i2][j1][j2]) ||
                    !isFinite(tstack2Enthalpies[i1][i2][j1][j2])    )
                {
                    tstack2Entropies [i1][i2][j1][j2] = -1.0;
                    tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
                }
            }
    }

    void getTriloop()
    {
        triloopEntropies .clear();  // ?? dont reuse
        triloopEnthalpies.clear();

        triloop_ds->seekg(std::ios_base::beg);
        triloop_dh->seekg(std::ios_base::beg);

        while ( readTLoop(triloop_ds, triloopEntropies ) ;
        while ( readTLoop(triloop_dh, triloopEnthalpies) ;
    }

    void getTetraloop()
    {
        tetraloopEntropies .clear();  // ?? dont reuse
        tetraloopEnthalpies.clear();

        triloop_ds->seekg(std::ios_base::beg);
        triloop_dh->seekg(std::ios_base::beg);

        while ( readTLoop(triloop_ds, tetraloopEntropies ) ;
        while ( readTLoop(triloop_dh, tetraloopEnthalpies) ;
    }

    void tableStartATS(double atp_value)
    {
        for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            atpS[i][j] = 0.00000000001;

        atpS[0][3] = atpS[3][0] = atp_value;
    }

    void tableStartATH(double atp_value)
    {
        for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            atpH[i][j] = 0.0;

        atpH[0][3] = atpH[3][0] = atp_value;
    }

};

thal_parameters::thal_parameters()
    : m_thal_p  { std::make_unique<thal_parameters::impl>() }
{}

thal_parameters::thal_parameters(const std::filesystem::path& dirname )
    : m_thal_p  { std::make_unique<thal_parameters::impl>(const std::filesystem::path& dirname ) }
{}

thal_parameters::~thal_parameters()
{}

int thal_parameters::set_defaults( )
{
    this->m_thal_p->set_defaults();
    return 0;
}

int thal_parameters::load (const std::filesystem::path& dirname )
{
    this->m_thal_p->load(dirname);
    return 0;
}






//                    ***************   CProgParam_ThAl   ***********************


/* Set default args for oligo */
void
CProgParam_ThAl::set_oligo_defaults( )

//                    ***************   ThAl impl  ***********************
#ifndef MIN_HRPN_LOOP
#define MIN_HRPN_LOOP 3 /*  minimum size of hairpin loop */
#endif

class ThAl
{
    static constexpr double  SMALL_NON_ZERO {0.000001};
    static int DBL_EQ(double X, double Y) { return (X - Y) < SMALL_NON_ZERO ? 1 : 2; /* 1 when numbers are equal */}
    static int max5(double a, double b, double c, double d, double e)
    {
        if(a > b &&
           a > c &&
           a > d &&
           a > e    ) return 1;
        else if(b > c &&
                b > d &&
                b > e    ) return 2;
        else if(c > d &&
                c > e    ) return 3;
        else if(d > e    ) return 4;
        else return 5;
    }

    static constexpr double R    = 1.9872;          /* cal/Kmol */
    static constexpr double ILAS = (-300 / 310.15); /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
    static constexpr double ILAH = 0.0;             /* Internal Loop EntHalpy Asymmetry correction */
    static constexpr double AT_H = 2200.0;          /* AT penalty */
    static constexpr double AT_S = 6.9;             /* AT penalty */
    static constexpr double MinEntropyCutoff = -2500.0; /* to filter out non-existing entropies */
    static constexpr double MinEntropy       = -3224.0; /* initiation */
    static constexpr double G2               = 0.0; /* structures w higher G are considered to be unstabile */
    static constexpr ABSOLUTE_ZERO    = 273.15;
    static constexpr int    MIN_LOOP         = 0;

    /* w/o init not constant anymore, cause for unimolecular and bimolecular foldings there are different values */
    double dplx_init_H;    /* initiation enthalpy; for duplex 200, for unimolecular structure 0 */
    double dplx_init_S;    /* initiation entropy; for duplex -5.7, for unimoleculat structure 0 */
    double saltCorrection; /* value calculated by saltCorrectS,
                                     includes correction for monovalent and divalent cations */
    double RC;             /* universal gas constant multiplied w DNA conc - for melting temperature */
    double SHleft;         /* var that helps to find str w highest melting temperature TODO not initialized !!!!??  */
    int bestI, bestJ;      /* starting position of most stable str */

    seq            oligo1,  oligo2;   /* inserted oligo sequenced */
    seq            numSeq1, numSeq2;  /* same as oligo1 and oligo2 but converted to numbers */
    static int    len1, len2, len3;   /* length of sequense 1 and 2 */

    struct tracer /* structure only for tracebacku - unimolecular str */
    {
        int i;
        int j;
        int mtrx; /* [0 1] EntropyDPT/EnthalpyDPT*/
        struct tracer* next;
    };
    static void push(struct tracer** stack, int i, int j, int mtrx)
    {
        struct tracer* new_top;
        new_top = (struct tracer*) safe_malloc(sizeof(struct tracer), o);
        new_top->i = i;
        new_top->j = j;
        new_top->mtrx = mtrx;
        new_top->next = *stack;
        *stack = new_top;
    }

    std::vector<double> enthalpyDPT, entropyDPT,    /* DyProg matrix for values of enthalpy and entropy */
                        send5,       hend5;         /* calc 5'  */

    /// table where bp-s enthalpies, that retrieve to the most stable Tm, are saved
    double & EnthalpyDPT(int i, int j)   {  return enthalpyDPT[(j) + ((i-1)*len3) - (1)];  }

    /// table where bp-s entropies, that retrieve to the most stable Tm, are saved
    double & EntropyDPT (int i, int j)   {  return entropyDPT [(j) + ((i-1)*len3) - (1)];  }

    /// entropies of most stable hairpin terminal bp
    double & SEND5 (int i)   {  return send5[i];  }

    /// enthalpies of most stable hairpin terminal bp
    double & HEND5 (int i)   {  return hend5[i];  }

    bool & bpIndx (unsigned char a, unsigned char b)   {  return BPI[a][b];  } /* for traceing matrix BPI ?? */
    double & atPenaltyS (unsigned char a, unsigned char b)   {  return atpS[a][b];  }
    double & atPenaltyH (unsigned char a, unsigned char b)   {  return atpH[a][b];  }

    void initMatrix() ///< initiates DPTables of entropy and enthalpy for dimer
    {
        for (int i = 1; i <= len1; ++i)                         // 1 !!
        for (int j = 1; j <= len2; ++j)                         // 1 !!
            if (bpIndx(numSeq1[i], numSeq2[j]) == 0)            // bpIndx(a, b) BPI[a][b]
            {              //  mismatch
                EnthalpyDPT(i, j) = _INFINITY;                  // enthalpyDPT[(j) + ((i-1)*len3) - (1)]
                EntropyDPT (i, j) = -1.0;
            } else
            {              // watson crick match
                EnthalpyDPT(i, j) = 0.0;
                EntropyDPT (i, j) = MinEntropy;
            }
    }

    void initMatrix2() ///< initiates DPTables of entropy and enthalpy for monomer
    {
        for (int i = 1; i <= len1; ++i)                       // 1 !!
        for (int j = i; j <= len2; ++j)                       // 1 !!
            if ( j - i < MIN_HRPN_LOOP + 1 ||                 // i, j too near, loop too short
                 (bpIndx(numSeq1[i], numSeq1[j]) == 0))        // bpIndx(a, b) BPI[a][b]
            {            //  mismatch
                EnthalpyDPT(i, j) = _INFINITY;                 // enthalpyDPT[(j) + ((i-1)*len3) - (1)]
                EntropyDPT (i, j) = -1.0;
            } else
            {              // watson crick match
                EnthalpyDPT(i, j) = 0.0;
                EntropyDPT (i, j) = MinEntropy;
            }
    }
        public:

};



static void maxTM(int i, int j); /* finds max Tm while filling the dyn progr table using stacking S and stacking H (dimer) */
static void maxTM2(int i, int j); /* finds max Tm while filling the dyn progr table using stacking S and stacking H (monomer) */

static double Ss(int i, int j, int k); /* returns stack entropy */
static double Hs(int i, int j, int k); /* returns stack enthalpy */

/* traceback for dimers */
static void traceback(int i, int j, double RT, int* ps1, int* ps2, int maxLoop, thal_results* o);

/* traceback for hairpins */
static void tracebacku(int*, int, thal_results*);

/* prints ascii output of dimer structure */
char *drawDimer(int*, int*, double, double, double, const thal_mode mode, double, thal_results *);

static void save_append_string(char** ret, int *space, thal_results *o, const char *str);

static void save_append_char(char** ret, int *space, thal_results *o, const char str);

static int equal(double a, double b);

static void strcatc(char*, char);

/* terminal bp for monomer structure */
static void calc_terminal_bp(double temp);

/* executed in calc_terminal_bp; to find structure that corresponds to max Tm for terminal bp */
static double END5_1(int,int); /* END5_1(X,1/2) - 1=Enthalpy, 2=Entropy*/
static double END5_2(int,int);
static double END5_3(int,int);
static double END5_4(int,int);

static double Hd5(int,int); /* returns thermodynamic value (H) for 5' dangling end */
static double Hd3(int,int); /* returns thermodynamic value (H) for 3' dangling end */
static double Sd5(int,int); /* returns thermodynamic value (S) for 5' dangling end */
static double Sd3(int,int); /* returns thermodynamic value (S) for 3' dangling end */
static double Ststack(int,int); /* returns entropy value for terminal stack */
static double Htstack(int,int); /* returns enthalpy value for terminal stack */




static void 
fillMatrix(int maxLoop)                   /* calc-s thermod values into dynamic progr table (dimer) */
{
   int d, i, j, ii, jj;
   double SH[2];

   for (i = 1; i <= len1; ++i) {                                  // 1 !! because numSeq1[i-1] will be used
      for (j = 1; j <= len2; ++j) {                               // 1 !! numSeq2[j-1] used
         if(isFinite(EnthalpyDPT(i, j))) {                        /* if finite */
            SH[0] = -1.0;
            SH[1] = _INFINITY;
            LSH(i,j,SH);   /* calculate terminal S and terminal H reading from 5'end (Left hand/3' end - Right end) */
            if(isFinite(SH[1])) {
               EntropyDPT (i,j) = SH[0];
               EnthalpyDPT(i,j) = SH[1];                        // enthalpyDPT[(j) + ((i-1)*len3) - (1)]
            }
            if (i > 1 && j > 1) {
               maxTM(i, j); /* stack: sets EntropyDPT(i, j) and EnthalpyDPT(i, j) */
               for(d = 3; d <= maxLoop + 2; d++) { /* max=30, length over 30 is not allowed */
                  ii = i - 1;
                  jj = - ii - d + (j + i);
                  if (jj < 1) {
                     ii -= abs(jj-1);
                     jj = 1;
                  }
                  for (; ii > 0 && jj < j; --ii, ++jj) {
                     if (isFinite(EnthalpyDPT(ii, jj))) {
                        SH[0] = -1.0;
                        SH[1] = _INFINITY;     /* calculates bulges and internal loops for dimer structures */
                        calc_bulge_internal(ii, jj, i, j, SH,0,maxLoop);
                        if(SH[0] < MinEntropyCutoff) {
                           /* to not give dH any value if dS is unreasonable */
                           SH[0] = MinEntropy;
                           SH[1] = 0.0;
                        }
                        if(isFinite(SH[1])) {
                           EnthalpyDPT(i, j) = SH[1];
                           EntropyDPT(i, j) = SH[0];
                        }
                     }
                  }
               }
            } /* if */
         }
      } /* for */
   } /* for */
}

static void 
fillMatrix2(int maxLoop )         ///<  calc-s thermod values into dynamic progr table (monomer)
{
   int i, j;
   double SH [2];
   for (j = 2; j <= len2; ++j)
     for (i = j - MIN_HRPN_LOOP - 1; i >= 1; --i)
     {
        if (isFinite(EnthalpyDPT(i, j)))
        {
           SH[0] = -1.0;
           SH[1] = _INFINITY;
           maxTM2(i,j);                       /* calculate stack */
           CBI(i, j, SH, 0,maxLoop);          /* calculate Bulge and Internal loop and stack */
           SH[0] = -1.0;
           SH[1] = _INFINITY;
           calc_hairpin(i, j, SH, 0);
           if(isFinite(SH[1])) {
              if(SH[0] < MinEntropyCutoff){ /* to not give dH any value if dS is unreasonable */
                 SH[0] = MinEntropy;
                 SH[1] = 0.0;
              }
              EntropyDPT (i,j) = SH[0];
              EnthalpyDPT(i,j) = SH[1];
           }
        }
     }
}


static void 
maxTM(int i, int j)
{
   double T0, T1;
   double S0, S1;
   double H0, H1;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double), 0);
   T0 = T1 = -_INFINITY;
   S0 = EntropyDPT(i, j);
   H0 = EnthalpyDPT(i, j);
   RSH(i,j,SH);
   T0 = (H0 + dplx_init_H + SH[1]) /(S0 + dplx_init_S + SH[0] + RC); /* at current position */
   if(isFinite(EnthalpyDPT(i - 1, j - 1)) && isFinite(Hs(i - 1, j - 1, 1))) {
      S1 = (EntropyDPT(i - 1, j - 1) + Ss(i - 1, j - 1, 1));
      H1 = (EnthalpyDPT(i - 1, j - 1) + Hs(i - 1, j - 1, 1));
      T1 = (H1 + dplx_init_H + SH[1]) /(S1 + dplx_init_S + SH[0] + RC);
   } else {
      S1 = -1.0;
      H1 = _INFINITY;
      T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
   }
   
   if(S1 < MinEntropyCutoff) {
      /* to not give dH any value if dS is unreasonable */
      S1 = MinEntropy;
      H1 = 0.0;
   }
   if(S0 < MinEntropyCutoff) {
      /* to not give dH any value if dS is unreasonable */
      S0 = MinEntropy;
      H0 = 0.0;
   }
   if(T1 > T0) { 
      EntropyDPT(i, j) = S1;
      EnthalpyDPT(i, j) = H1;
   } else if(T0 >= T1) {
      EntropyDPT(i, j) = S0;
      EnthalpyDPT(i, j) = H0;
   }
   free(SH);
}

static void 
maxTM2(int i, int j) /* finds max Tm while filling the dyn progr table using stacking S and stacking H (monomer) */
{
   double T0, T1;
   double S0, S1;
   double H0, H1;
   T0 = T1 = -_INFINITY;
   S0 = EntropyDPT(i, j);
   H0 = EnthalpyDPT(i, j);
   T0 = (H0 + dplx_init_H) /(S0 + dplx_init_S + RC);    // ??
   if(isFinite(EnthalpyDPT(i, j)))
   {
      S1 = (EntropyDPT (i + 1, j - 1) + Ss(i, j, 2));
      H1 = (EnthalpyDPT(i + 1, j - 1) + Hs(i, j, 2));
   } else {
      S1 = -1.0;
      H1 = _INFINITY;
   }
   T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
   if(S1 < MinEntropyCutoff)
   {
      S1 = MinEntropy;
      H1 = 0.0;
   }
   if(S0 < MinEntropyCutoff)
   {
      S0 = MinEntropy;
      H0 = 0.0;
   }

   if(T1 > T0)
   {
      EntropyDPT (i, j) = S1;
      EnthalpyDPT(i, j) = H1;
   }
   else
   {
      EntropyDPT (i, j) = S0;
      EnthalpyDPT(i, j) = H0;
   }
}

/* calculate terminal entropy S and terminal enthalpy H starting reading from 5'end (Left hand/3' end - Right end) */
static void 
LSH(int i, int j, double EntropyEnthalpy[2])
{
   double S1, H1, T1, G1;
   double S2, H2, T2, G2;
   S1 = S2 = -1.0;
   H1 = H2 = -_INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {      // bpIndx(a, b) BPI[a][b] , mismatch
      EntropyDPT(i, j) = -1.0;
      EnthalpyDPT(i, j) = _INFINITY;               // enthalpyDPT[(j) + ((i-1)*len3) - (1)]
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(!isFinite(H1) || G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }
   /** If there is two dangling ends at the same end of duplex **/
   if((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1 ) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
        dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
        dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1<T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0) {
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   } else if ((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1<T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if (G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   } else if ((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;     
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2  && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0) {
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   G1 = H1 -TEMP_KELVIN*S1;   
   G2 = H2 -TEMP_KELVIN*S2;
   if(isFinite(H1)) {
      if(T1 < T2) {
         EntropyEnthalpy[0] = S2;
         EntropyEnthalpy[1] = H2;
      } else {
         EntropyEnthalpy[0] = S1;
         EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}
/* calculate terminal entropy S and terminal enthalpy H starting reading from 5'end (Left hand/3' end - Right end) */
static void
RSH(int i, int j, double* EntropyEnthalpy)
{
   double G1, G2;
   double S1, S2;
   double H1, H2;
   double T1, T2;
   S1 = S2 = -1.0;
   H1 = H2 = _INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(!isFinite(H1) || G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }
   
   if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
    G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }

   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2 >0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }

   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if (G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   G1 = H1 -TEMP_KELVIN*S1;
   G2 =  H2 -TEMP_KELVIN*S2;
   if(isFinite(H1)) {
      if(T1 < T2) {
         EntropyEnthalpy[0] = S2;
         EntropyEnthalpy[1] = H2;
      } else {
         EntropyEnthalpy[0] = S1;
         EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}

static double 
Ss(int i, int j, int k)
{
   if(k==2) {
      if (i >= j)
        return -1.0;
      if (i == len1 || j == len2 + 1)
        return -1.0;

      if (i > len1)
        i -= len1;
      if (j > len2)
        j -= len2;
      return stackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
   } else {
      return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}


static double 
Hs(int i, int j, int k)
{
   if(k==2) {
      if (i >= j)
        return _INFINITY;
      if (i == len1 || j == len2 + 1)
        return _INFINITY;

      if (i > len1)
        i -= len1;
      if (j > len2)
        j -= len2;
      if(isFinite(stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]])) {
         return stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
      } else {
         return _INFINITY;
      }
   } else {
      return stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}
/* carries out Bulge and Internal loop and stack calculations to hairpin */
static void
CBI(int i, int j, double EntropyEnthalpy[2], int traceback, int maxLoop)
{
   int d, ii, jj;
   for (d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop; --d)  // ??
     for (ii = i + 1; ii < j - d && ii <= len1; ++ii) {
        jj = d + ii;
        if(traceback==0)
        {
           EntropyEnthalpy[0] = -1.0;
           EntropyEnthalpy[1] = _INFINITY;
        }
        if (isFinite(EnthalpyDPT(ii, jj)) && isFinite(EnthalpyDPT(i, j)))
        {
           calc_bulge_internal2(i, j, ii, jj, EntropyEnthalpy, traceback,maxLoop);
           if(isFinite(EntropyEnthalpy[1])) {
              if(EntropyEnthalpy[0] < MinEntropyCutoff) {
                 EntropyEnthalpy[0] = MinEntropy;
                 EntropyEnthalpy[1] = 0.0;
              }
              if(traceback==0) {
                 EnthalpyDPT(i, j) = EntropyEnthalpy[1];
                 EntropyDPT(i, j) = EntropyEnthalpy[0];
              }
           }
        }
     }
   return;
}
/* finds monomer structure that has maximum Tm */
static void
calc_hairpin(int i, int j, double EntropyEnthalpy[2], int traceback)
{
   int loopSize = j - i - 1;
   double G1, G2;
   G1 = G2 = -_INFINITY;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double), 0);
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   if(loopSize < MIN_HRPN_LOOP) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   if (i <= len1 && len2 < j) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   } else if (i > len2) {
      i -= len1;
      j -= len2;
   }
   if(loopSize <= 30) {
      EntropyEnthalpy[1] = hairpinLoopEnthalpies[loopSize - 1];
      EntropyEnthalpy[0] = hairpinLoopEntropies[loopSize - 1];
   } else {
      EntropyEnthalpy[1] = hairpinLoopEnthalpies[29];
      EntropyEnthalpy[0] = hairpinLoopEntropies[29];
   }

   if (loopSize > 3) { /* for loops 4 bp and more in length, terminal mm are accounted */
      EntropyEnthalpy[1] += tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
      EntropyEnthalpy[0] += tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
   } else if(loopSize == 3){ /* for loops 3 bp in length at-penalty is considered */
      EntropyEnthalpy[1] += atPenaltyH(numSeq1[i], numSeq1[j]);
      EntropyEnthalpy[0] += atPenaltyS(numSeq1[i], numSeq1[j]);
   }

   if (loopSize == 3) {         /* closing AT-penalty (+), triloop bonus, hairpin of 3 (+) */
      struct triloop* loop;
      if (numTriloops)
      {
         if ((loop = (struct triloop*) bsearch( numSeq1 + i,
                                                triloopEnthalpies,
                                                numTriloops,
                                                sizeof(struct triloop),
                                                comp3loop))                   )
           EntropyEnthalpy[1] += loop->value;

         if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEntropies, numTriloops, sizeof(struct triloop), comp3loop)))
           EntropyEnthalpy[0] += loop->value;
      }
   } else if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
      struct tetraloop* loop;
      if (numTetraloops) {
         if ((loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEnthalpies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
            EntropyEnthalpy[1] += loop->value;
         }
         if ((loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEntropies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
            EntropyEnthalpy[0] += loop->value;
         }
      }
   }
   if(!isFinite(EntropyEnthalpy[1])) {
      EntropyEnthalpy[1] = _INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   if(isPositive(EntropyEnthalpy[1]) && isPositive(EntropyEnthalpy[0]) && (!isPositive(EnthalpyDPT(i, j)) || !isPositive(EntropyDPT(i, j)))) { /* if both, S and H are positive */
      EntropyEnthalpy[1] = _INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   RSH(i,j,SH);
   G1 = EntropyEnthalpy[1]+SH[1] -TEMP_KELVIN*(EntropyEnthalpy[0]+SH[0]);
   G2 = EnthalpyDPT(i, j)+SH[1] -TEMP_KELVIN*(EntropyDPT(i, j)+SH[0]);
     if(G2 < G1 && traceback == 0) {
      EntropyEnthalpy[0] = EntropyDPT(i, j);
      EntropyEnthalpy[1] = EnthalpyDPT(i, j);
   }
   free(SH);
   return;
}

/* calculates bulges and internal loops for dimer structures */
static void 
calc_bulge_internal(int i, int j, int ii, int jj, double EntropyEnthalpy[2], int traceback, int maxLoop)
{
   int loopSize1, loopSize2, loopSize;
   double S,H,G1,G2;
   int N, N_loop;
   double SH[2];
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   S = -1.0;
   H = _INFINITY;
   loopSize1 = ii - i - 1;
   loopSize2 = jj - j - 1;
   if(ii < jj) {
      N = ((2 * i)/2);
      N_loop = N;
      if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
      if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
   } else {
      N = ((2 * j)/2);
      N_loop = 2 * jj;
      if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
      if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
      N_loop = (N_loop/2) - 1;
   }
#ifdef DEBUG
   if (ii <= i){
      fputs("Error in calc_bulge_internal(): ii is not greater than i\n", stderr);
   }
   if (jj <= j)
     fputs("Error in calc_bulge_internal(): jj is not greater than j\n", stderr);
#endif

#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      free(SH);
      return;
   }
#endif
#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      free(SH);
      return;
   }
#endif
   loopSize = loopSize1 + loopSize2-1;
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                              the intervening nn-pair must be added */

         if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
            H = bulgeLoopEnthalpies[loopSize] +
              stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
            S = bulgeLoopEntropies[loopSize] +
              stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
         }
         if(isPositive(H) || isPositive(S)){
            H = _INFINITY;
            S = -1.0;
         }
         H += EnthalpyDPT(i, j);
         S += EntropyDPT(i, j);
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
         RSH(ii,jj,SH);
         G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
         G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*((EntropyDPT(ii, jj)+SH[0]));
         if((G1< G2) || (traceback==1)) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

         H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
         H += EnthalpyDPT(i, j);

         S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
         S += EntropyDPT(i, j);
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
         if(isPositive(H) && isPositive(S)){ 
            H = _INFINITY;
            S = -1.0;
         }
        
     RSH(ii,jj,SH);
         G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
         G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
         if(G1< G2 || (traceback==1)){
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
         
      }
   } else if (loopSize1 == 1 && loopSize2 == 1) {
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        stackint2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      S += EntropyDPT(i, j);

      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        stackint2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      H += EnthalpyDPT(i, j);
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
      if(isPositive(H) && isPositive(S)){
         H = _INFINITY;
         S = -1.0;
      }    
     RSH(ii,jj,SH);
         G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
      G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
           if((G1< G2) || traceback==1) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      free(SH);
      return;
   } else { /* only internal loops */
      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        tstackEnthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]
        + (ILAH * abs(loopSize1 - loopSize2));
      H += EnthalpyDPT(i, j);

      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        tstackEntropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      S += EntropyDPT(i, j);
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
   if(isPositive(H) && isPositive(S)){ 
         H = _INFINITY;
         S = -1.0;
      }
     RSH(ii,jj,SH);
     G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
     G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
     if((G1< G2) || (traceback==1)){
             EntropyEnthalpy[0] = S;
             EntropyEnthalpy[1] = H;
      }
   }
   return;
}
/* calculates bulges and internal loops for monomer structures */
static void 
calc_bulge_internal2(int i, int j, int ii, int jj, double EntropyEnthalpy[2], int traceback, int maxLoop)
{
   int loopSize1, loopSize2, loopSize;
   double T1, T2;
   double S,H;
   /* int N, N_loop; Triinu, please review */
   T1 = T2 = -_INFINITY;
   S = MinEntropy;
   H = 0.0;
   loopSize1 = ii - i - 1;
   loopSize2 = j - jj - 1;
   if (loopSize1 + loopSize2 > maxLoop) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   /* Triinu, please review the statements below. */
   /* if(i < (len1 -j)) { */
     /* N  = i; */
      /* N_loop = (i - 1); */
   /* } else { */
     /* N = len1-j;  */
      /* N_loop = len1 - j - 1; */
   /* } */
#ifdef DEBUG
   if (ii <= i)
     fputs("Error in calc_bulge_internal(): ii isn't greater than i\n", stderr);
   if (jj >= j)
     fputs("Error in calc_bulge_internal(): jj isn't less than j\n", stderr);
   if (ii >= jj)
     fputs("Error in calc_bulge_internal(): jj isn't greater than ii\n", stderr);

   if ((i <= len1 && len1 < ii) || (jj <= len2 && len2 < j))  {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
#endif

#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      return;
   }
#endif
#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      return;
   }
#endif

#ifdef DEBUG
   if (i > len1)
     i -= len1;
   if (ii > len1)
     ii -= len1;
   if (j > len2)
     j -= len2;
   if (jj > len2)
     jj -= len2;
#endif
   loopSize = loopSize1 + loopSize2 -1; /* for indx only */
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                              the intervening nn-pair must be added */
         if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
            H = bulgeLoopEnthalpies[loopSize] +
              stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
            S = bulgeLoopEntropies[loopSize] +
              stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
         }
         if(traceback!=1) {
            H += EnthalpyDPT(ii, jj); /* bulge koos otsaga, st bulge i,j-ni */
            S += EntropyDPT(ii, jj);
         }

         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
        
         T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
         T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);

         if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }

      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

         H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
         if(traceback!=1)                      //  atpH[a][b]
           H += EnthalpyDPT(ii, jj);

         S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
         if(traceback!=1)
           S += EntropyDPT(ii, jj);
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
          
         T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
         T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);

         if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      }
   } /* end of calculating bulges */
   else if (loopSize1 == 1 && loopSize2 == 1) {
      /* mismatch nearest neighbor parameters */

      S = stackint2Entropies[numSeq1[i]] [numSeq1[i+1 ]][numSeq2[j ]][numSeq2[j-1]] +
          stackint2Entropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
        S += EntropyDPT(ii, jj);

      H = stackint2Enthalpies[numSeq1[i ]][numSeq1[i+1 ]][numSeq2[j ]][numSeq2[j-1 ]] +
          stackint2Enthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
        H += EnthalpyDPT(ii, jj);
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
    
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);

      if((DBL_EQ(T1,T2) == 2) || traceback) {
         if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      }
      return;
   } else { /* only internal loops */

      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i ]][numSeq1[i+1 ]][numSeq2[ j]][numSeq2[j-1 ]] +
                                             tstackEnthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]]
        + (ILAH * abs(loopSize1 - loopSize2));
      if(traceback!=1)
        H += EnthalpyDPT(ii, jj);

      S = interiorLoopEntropies[loopSize] +
              tstackEntropies[numSeq1[i ]][numSeq1[i+1 ]][numSeq2[j ]][numSeq2[j-1 ]] +
              tstackEntropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      if(traceback!=1)
        S += EntropyDPT(ii, jj);
      if(!isFinite(H))
      {
         H = _INFINITY;
         S = -1.0;
      }
    
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);
      if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
         EntropyEnthalpy[0] = S;
         EntropyEnthalpy[1] = H;
      }
   }
   return;
}

static void 
calc_terminal_bp(double temp) /* compute exterior loop */ /* terminal bp for monomer structure */
{
   int i;
   int max;
   SEND5(0) = SEND5(1) = -1.0;
   HEND5(0) = HEND5(1) = _INFINITY;
   for(i = 2; i<=(len1); i++) {
      SEND5(i) = MinEntropy;
      HEND5(i) = 0;
   }

   double T1, T2, T3, T4, T5;
   T1 = T2 = T3 = T4 = T5 = -_INFINITY;
   double G;
   /* adding terminal penalties to 3' end and to 5' end */
   for(i = 2; i <= len1; ++i) {
      max = 0;
      T1 = T2 = T3 = T4 = T5 = -_INFINITY;
      T1 = (HEND5 (i-1) + dplx_init_H) / (SEND5(i- 1) + dplx_init_S + RC);
      T2 = (END5_1(i,1) + dplx_init_H) / (END5_1(i,2) + dplx_init_S + RC);
      T3 = (END5_2(i,1) + dplx_init_H) / (END5_2(i,2) + dplx_init_S + RC);
      T4 = (END5_3(i,1) + dplx_init_H) / (END5_3(i,2) + dplx_init_S + RC);
      T5 = (END5_4(i,1) + dplx_init_H) / (END5_4(i,2) + dplx_init_S + RC);
      max = max5(T1,T2,T3,T4,T5);
      switch (max) {
       case 1:
         SEND5(i) = SEND5(i - 1);
         HEND5(i) = HEND5(i - 1);
         break;
       case 2:
         G = END5_1(i,1) - (temp * (END5_1(i,2)));
         if(G < G2) {
            SEND5(i) = END5_1(i,2);
            HEND5(i) = END5_1(i,1);
         } else {
            SEND5(i) = SEND5(i - 1);
            HEND5(i) = HEND5(i - 1);
         }
         break;
       case 3:
         G = END5_2(i,1) - (temp * (END5_2(i,2)));
         if(G < G2) {
            SEND5(i) = END5_2(i,2);
            HEND5(i) = END5_2(i,1);
         } else {
            SEND5(i) = SEND5(i - 1);
            HEND5(i) = HEND5(i - 1);
         }
         break;
       case 4:
         G = END5_3(i,1) - (temp * (END5_3(i,2)));
         if(G < G2) {
            SEND5(i) = END5_3(i,2);
            HEND5(i) = END5_3(i,1);
         } else {
            SEND5(i) = SEND5(i - 1);
            HEND5(i) = HEND5(i - 1);
         }
         break;
       case 5:
         G = END5_4(i,1) - (temp * (END5_4(i,2)));
         if(G < G2) {
            SEND5(i) = END5_4(i,2);
            HEND5(i) = END5_4(i,1);
         } else {
            SEND5(i) = SEND5(i - 1);
            HEND5(i) = HEND5(i - 1);
         }
         break;
       default:
#ifdef DEBUG
         printf ("WARNING: max5 returned character code %d ??\n", max);
#endif
         break;
      }
   }
}

static double 
END5_1(int i,int hs)
{
   int k;
   double max_tm; /* energy min */
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   S_max = S = -1.0;
   T1 = T2 = -_INFINITY;
   max_tm = -_INFINITY;
   for(k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
         H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
         S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
         if(!isFinite(H) || H > 0 || S > 0) { /* H and S must be greater than 0 to avoid BS */
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
         H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
         S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
         if(S > MinEntropyCutoff) {
            H_max = H;
            S_max = S;
            max_tm = T1;
         }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double 
END5_2(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
         H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
         S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
         H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
         S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
         if(S > MinEntropyCutoff) {
            H_max = H;
            S_max = S;
            max_tm = T1;
         }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double 
END5_3(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
         H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
         S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
         H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
         S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
         if(S > MinEntropyCutoff) {
            H_max = H;
            S_max = S;
            max_tm = T1;
         }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double 
END5_4(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for(k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
         H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
         S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
         H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
         S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
         if(!isFinite(H) || H > 0 || S > 0) {
            H = _INFINITY;
            S = -1.0;
         }
         T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
         if(S > MinEntropyCutoff) {
            H_max = H;
            S_max = S;
            max_tm = T1;
         }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}


static double 
Sd5(int i, int j)
{
   return dangleEntropies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
}

static double 
Hd5(int i, int j)
{
   return dangleEnthalpies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
}

static double 
Sd3(int i, int j)
{
   return dangleEntropies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
}

static double 
Hd3(int i, int j)
{
   return dangleEnthalpies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
}

static double 
Ststack(int i, int j)
{
   return tstack2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
}

static double 
Htstack(int i, int j)
{ /* e.g AG_TC 210 */
   return tstack2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
}
/* Is sequence symmetrical */
/* Return if string is symmetrical. Used only once only in thal. TODO change to use code */
static bool
symmetry_thermo(const seq& sq)
{
   int end= sq.length();
   if(end%2)  return false;   // symmetrical only if 2 | L
   end--;
   for(int i=0; i<end; i++, end-- )
   {
      char s=nt2code(sq[i]  );                  //  to avoid this
      char e=nt2code(sq[end]);
      if (! bpIndx(s,e) )  return false;
   }
   return true;
}

static void 
tracebacku(std::vector<int>& bp, int maxLoop,thal_results* o) /* traceback for unimolecular structure */ /* traceback for hairpins */
{
   int i, j;
   i = j = 0;
   int ii, jj, k;
   struct tracer *top, *stack = NULL;
   double* SH1;
   double* SH2;
   double* EntropyEnthalpy;
   SH1 = (double*) safe_malloc(2 * sizeof(double), o);
   SH2 = (double*) safe_malloc(2 * sizeof(double), o);
   EntropyEnthalpy = (double*) safe_malloc(2 * sizeof(double), o);
   push(&stack,len1, 0, 1, o);
   while(stack) {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;
      if(top->mtrx==1) {
         while (equal(SEND5(i), SEND5(i - 1)) && equal(HEND5(i), HEND5(i - 1))) /* if previous structure is the same as this one */
           --i;
         if (i == 0)
           continue;
         if (equal(SEND5(i), END5_1(i,2)) && equal(HEND5(i), END5_1(i,1))) {
            for (k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k)
              if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i)) &&
                  equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i))) {
                 push(&stack, k + 1, i,0, o);
                 break;
              }
            else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i)) &&
                     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i))) {
               push(&stack, k + 1, i, 0, o);
               push(&stack, k, 0, 1, o);
               break;
            }
         }
         else if (equal(SEND5(i), END5_2(i,2)) && equal(HEND5(i), END5_2(i,1))) {
            for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
              if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i)) &&
                  equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i))) {
                 push(&stack, k + 2, i, 0, o);
                 break;
              }
            else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i)) &&
                     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i))) {
               push(&stack, k + 2, i, 0, o);
               push(&stack, k, 0, 1, o);
               break;
            }
         }
         else if (equal(SEND5(i), END5_3(i,2)) && equal(HEND5(i), END5_3(i,1))) {
            for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
              if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1))
                  && equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
                 push(&stack, k + 1, i - 1, 0, o);
                 break;
              }
            else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1)) &&
                     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
               push(&stack, k + 1, i - 1, 0, o); /* matrix 0  */
               push(&stack, k, 0, 1, o); /* matrix 3 */
               break;
            }
         }
         else if(equal(SEND5(i), END5_4(i,2)) && equal(HEND5(i), END5_4(i,1))) {
            for (k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k)
              if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
                  equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1))) {
                 push(&stack, k + 2, i - 1, 0, o);
                 break;
              }
            else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
                     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1)) ) {
               push(&stack, k + 2, i - 1, 0, o);
               push(&stack, k, 0, 1, o);
               break;
            }
         }
      }
      else if(top->mtrx==0) {
         bp[i - 1] = j;
         bp[j - 1] = i;
         SH1[0] = -1.0;
         SH1[1] = _INFINITY;
         calc_hairpin(i, j, SH1, 1); /* 1 means that we use this method in traceback */
         SH2[0] = -1.0;
         SH2[1] = _INFINITY;
         CBI(i,j,SH2,2,maxLoop);
         if (equal(EntropyDPT(i, j), Ss(i, j, 2) + EntropyDPT(i + 1, j - 1)) &&
             equal(EnthalpyDPT(i, j), Hs(i, j, 2) + EnthalpyDPT(i + 1, j - 1))) {
            push(&stack, i + 1, j - 1, 0, o);
         }
         else if (equal(EntropyDPT(i, j), SH1[0]) && equal(EnthalpyDPT(i,j), SH1[1]));
         else if (equal(EntropyDPT(i, j), SH2[0]) && equal(EnthalpyDPT(i, j), SH2[1])) {
            int d, done;
            for (done = 0, d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop && !done; --d)
              for (ii = i + 1; ii < j - d; ++ii) {
                 jj = d + ii;
                 EntropyEnthalpy[0] = -1.0;
                 EntropyEnthalpy[1] = _INFINITY;
                 calc_bulge_internal2(i, j, ii, jj,EntropyEnthalpy,1,maxLoop);
                 if (equal(EntropyDPT(i, j), EntropyEnthalpy[0] + EntropyDPT(ii, jj)) &&
                     equal(EnthalpyDPT(i, j), EntropyEnthalpy[1] + EnthalpyDPT(ii, jj))) {
                    push(&stack, ii, jj, 0, o);
                    ++done;
                    break;
                 }
              }
         } else {
         }
      }
      free(top);
   }
   free(SH1);
   free(SH2);
   free(EntropyEnthalpy);
}

/* traceback for dimers */
static void 
traceback(int i, int j, double RT, std::vector<int>&  ps1, std::vector<int>& ps2, int maxLoop, thal_results* o)
{
   int d, ii, jj, done;
   double SH[2];
   ps1[i - 1] = j;
   ps2[j - 1] = i;
   while(1)
   {
      SH[0] = -1.0;
      SH[1] = _INFINITY;
      LSH(i,j,SH);
      if(equal(EntropyDPT(i,j),SH[0]) && equal(EnthalpyDPT(i,j),SH[1])) {
         break;
      }
      done = 0;
      if (i > 1 && j > 1 && equal(EntropyDPT (i,j), Ss(i - 1, j - 1, 1) + EntropyDPT (i - 1, j - 1)) &&
                            equal(EnthalpyDPT(i,j), Hs(i - 1, j - 1, 1) + EnthalpyDPT(i - 1, j - 1)))
      {
         i = i - 1;
         j = j - 1;
         ps1[i - 1] = j;
         ps2[j - 1] = i;
         done = 1;
      }
      for (d = 3; !done && d <= maxLoop + 2; ++d) {
         ii = i - 1;
         jj = -ii - d + (j + i);
         if (jj < 1) {
            ii -= abs(jj-1);
            jj = 1;
         }
         for (; !done && ii > 0 && jj < j; --ii, ++jj) {
            SH[0] = -1.0;
            SH[1] = _INFINITY;
            calc_bulge_internal(ii, jj, i, j, SH,1,maxLoop);
            if (equal(EntropyDPT (i, j), SH[0]) &&
                equal(EnthalpyDPT(i, j), SH[1]))
            {
               i = ii;
               j = jj;
               ps1[i - 1] = j;
               ps2[j - 1] = i;
               done = 1;
               break;
            }
         }
      }
   }
}

char * 
drawDimer(std::vector<int>& ps1, std::vector<int>& ps2, double temp, double H, double S, const CProgParam_ThAl::mode mode, double t37, thal_results *o)
{
   int  ret_space = 0;
   int ret_nr, ret_pr_once;
   std::string ret_str[4];
   int i, j, k, numSS1, numSS2, N;
   char* duplex[4];
   double G, t;
   t = G = 0;
   if (!isFinite(temp)){
      if((mode != THL_FAST) && (mode != THL_DEBUG_F) && (mode != THL_STRUCT)) {
         printf("No predicted secondary structures for given sequences\n");
      }
      o->temp = 0.0; /* lets use generalization here; this should rather be very negative value */
      strcpy(o->msg, "No predicted sec struc for given seq");
      return NULL;
   } else {
      N=0;
      for(i=0;i<len1;i++){
         if(ps1[i]>0) ++N;
      }
      for(i=0;i<len2;i++) {
         if(ps2[i]>0) ++N;
      }
      N = (N/2) -1;
      t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
         G = (H) - (t37 * (S + (N * saltCorrection)));
         S = S + (N * saltCorrection);
         o->temp = (double) t;
         /* maybe user does not need as precise as that */
         /* printf("Thermodynamical values:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\tN = %d, SaltC=%f, RC=%f\n",
                len1, (double) S, (double) H, (double) G, (double) t, (int) N, saltCorrection, RC); */
         if (mode != THL_STRUCT) {
           printf("Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
                  (double) S, (double) H, (double) G, (double) t);
         } else {
           sprintf(ret_para, "t: %.1f  dG: %.0f  dH: %.0f  dS: %.0f\\n",
                   (double) t, (double) G, (double) H, (double) S);
         }
      } else {
         o->temp = (double) t;
         return NULL;
      }
   }

   duplex[0] = (char*) safe_malloc(len1 + len2 + 1, o);
   duplex[1] = (char*) safe_malloc(len1 + len2 + 1, o);
   duplex[2] = (char*) safe_malloc(len1 + len2 + 1, o);
   duplex[3] = (char*) safe_malloc(len1 + len2 + 1, o);
   duplex[0][0] = duplex[1][0] = duplex[2][0] = duplex[3][0] = 0;

   i = 0;
   numSS1 = 0;
   while (ps1[i++] == 0) ++numSS1;
   j = 0;
   numSS2 = 0;
   while (ps2[j++] == 0) ++numSS2;

   if (numSS1 >= numSS2){
      for (i = 0; i < numSS1; ++i) {
         strcatc(duplex[0], oligo1[i]);
         strcatc(duplex[1], ' ');
         strcatc(duplex[2], ' ');
      }
      for (j = 0; j < numSS1 - numSS2; ++j) strcatc(duplex[3], ' ');
      for (j = 0; j < numSS2; ++j) strcatc(duplex[3], oligo2[j]);
   } else {
      for (j = 0; j < numSS2; ++j) {
         strcatc(duplex[3], oligo2[j]);
         strcatc(duplex[1], ' ');
         strcatc(duplex[2], ' ');
      }
      for (i = 0; i < numSS2 - numSS1; ++i)
        strcatc(duplex[0], ' ');
      for (i = 0; i < numSS1; ++i)
        strcatc(duplex[0], oligo1[i]);
   }
   i = numSS1 + 1;
   j = numSS2 + 1;

   while (i <= len1) {
      while (i <= len1 && ps1[i - 1] != 0 && j <= len2 && ps2[j - 1] != 0) {
         strcatc(duplex[0], ' ');
         strcatc(duplex[1], oligo1[i - 1]);
         strcatc(duplex[2], oligo2[j - 1]);
         strcatc(duplex[3], ' ');
         ++i;
         ++j;
      }
      numSS1 = 0;
      while (i <= len1 && ps1[i - 1] == 0) {
         strcatc(duplex[0], oligo1[i - 1]);
         strcatc(duplex[1], ' ');
         ++numSS1;
         ++i;
      }
      numSS2 = 0;
      while (j <= len2 && ps2[j - 1] == 0) {
         strcatc(duplex[2], ' ');
         strcatc(duplex[3], oligo2[j - 1]);
         ++numSS2;
         ++j;
      }
      if (numSS1 < numSS2)
        for (k = 0; k < numSS2 - numSS1; ++k) {
           strcatc(duplex[0], '-');
           strcatc(duplex[1], ' ');
        }
      else if (numSS1 > numSS2)
        for (k = 0; k < numSS1 - numSS2; ++k) {
           strcatc(duplex[2], ' ');
           strcatc(duplex[3], '-');
        }
   }
   if ((mode == THL_GENERAL) || (mode == THL_DEBUG)) {
     printf("SEQ\t");
     printf("%s\n", duplex[0]);
     printf("SEQ\t");
     printf("%s\n", duplex[1]);
     printf("STR\t");
     printf("%s\n", duplex[2]);
     printf("STR\t");
     printf("%s\n", duplex[3]);
   }
   if (mode == THL_STRUCT) {
     ret_str[3] = NULL;
     ret_str[0] = (char*) safe_malloc(len1 + len2 + 10, o);
     ret_str[1] = (char*) safe_malloc(len1 + len2 + 10, o);
     ret_str[2] = (char*) safe_malloc(len1 + len2 + 10, o);
     ret_str[0][0] = ret_str[1][0] = ret_str[2][0] = '\0';

     /* Join top primer */
     strcpy(ret_str[0], "   ");
     strcat(ret_str[0], duplex[0]);
     ret_nr = 0;
     while (duplex[1][ret_nr] != '\0') {
       if (duplex[1][ret_nr] == 'A' || duplex[1][ret_nr] == 'T' || 
           duplex[1][ret_nr] == 'C' || duplex[1][ret_nr] == 'G' || 
           duplex[1][ret_nr] == '-') {
         ret_str[0][ret_nr + 3] = duplex[1][ret_nr];
       }
       ret_nr++;
     }
     if (strlen(duplex[1]) > strlen(duplex[0])) {
       ret_str[0][strlen(duplex[1]) + 3] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[0]) - 1;
     while (ret_nr > 0 && (ret_str[0][ret_nr] == ' ' || ret_str[0][ret_nr] == '-')) {
       ret_str[0][ret_nr--] = '\0';
     }
     /* Write the 5' */
     ret_nr = 3;
     ret_pr_once = 1;
     while (ret_str[0][ret_nr] != '\0' && ret_pr_once == 1) {
       if (ret_str[0][ret_nr] == 'A' || ret_str[0][ret_nr] == 'T' ||
           ret_str[0][ret_nr] == 'C' || ret_str[0][ret_nr] == 'G' ||
           ret_str[0][ret_nr] == '-') {
         ret_str[0][ret_nr - 3] = '5';
	 ret_str[0][ret_nr - 2] = '\'';
	 ret_pr_once = 0;
       }
       ret_nr++;
     }

     /* Create the align tics */
     strcpy(ret_str[1], "     ");
     for (i = 0 ; i < strlen(duplex[1]) ; i++) {
       if (duplex[1][i] == 'A' || duplex[1][i] == 'T' || 
           duplex[1][i] == 'C' || duplex[1][i] == 'G' ) {
	 ret_str[1][i + 3] = '|';
       } else {
         ret_str[1][i + 3] = ' ';
       }
       ret_str[1][i + 4] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[1]) - 1;
     while (ret_nr > 0 && ret_str[1][ret_nr] == ' ') {
       ret_str[1][ret_nr--] = '\0';
     }
     /* Join bottom primer */
     strcpy(ret_str[2], "   ");
     strcat(ret_str[2], duplex[2]);
     ret_nr = 0;
     while (duplex[3][ret_nr] != '\0') {
       if (duplex[3][ret_nr] == 'A' || duplex[3][ret_nr] == 'T' ||
           duplex[3][ret_nr] == 'C' || duplex[3][ret_nr] == 'G' ||
           duplex[3][ret_nr] == '-') {
         ret_str[2][ret_nr + 3] = duplex[3][ret_nr];
       }
       ret_nr++;
     }
     if (strlen(duplex[3]) > strlen(duplex[2])) {
       ret_str[2][strlen(duplex[3]) + 3] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[2]) - 1;
     while (ret_nr > 0 && (ret_str[2][ret_nr] == ' ' || ret_str[2][ret_nr] == '-')) {
       ret_str[2][ret_nr--] = '\0';
     }
     /* Write the 5' */
     ret_nr = 3;
     ret_pr_once = 1;
     while (ret_str[2][ret_nr] != '\0' && ret_pr_once == 1) {
       if (ret_str[2][ret_nr] == 'A' || ret_str[2][ret_nr] == 'T' ||
           ret_str[2][ret_nr] == 'C' || ret_str[2][ret_nr] == 'G' ||
           ret_str[2][ret_nr] == '-') {
         ret_str[2][ret_nr - 3] = '3';
         ret_str[2][ret_nr - 2] = '\'';
         ret_pr_once = 0;
       }
       ret_nr++;
     }

     save_append_string(&ret_str[3], &ret_space, o, ret_para);
     save_append_string(&ret_str[3], &ret_space, o, ret_str[0]);
     save_append_string(&ret_str[3], &ret_space, o, " 3\'\\n");
     save_append_string(&ret_str[3], &ret_space, o, ret_str[1]);
     save_append_string(&ret_str[3], &ret_space, o, "\\n");
     save_append_string(&ret_str[3], &ret_space, o, ret_str[2]);
     save_append_string(&ret_str[3], &ret_space, o, " 5\'\\n");


/*
     save_append_string(&ret_str, &ret_space, o, "SEQ ");
     save_append_string(&ret_str, &ret_space, o, duplex[0]);
     save_append_string(&ret_str, &ret_space, o, "\\nSEQ ");
     save_append_string(&ret_str, &ret_space, o, duplex[1]);
     save_append_string(&ret_str, &ret_space, o, "\\nSTR ");
     save_append_string(&ret_str, &ret_space, o, duplex[2]);
     save_append_string(&ret_str, &ret_space, o, "\\nSTR ");
     save_append_string(&ret_str, &ret_space, o, duplex[3]);
     save_append_string(&ret_str, &ret_space, o, "\\n");
*/
     ret_ptr = (char *) safe_malloc(strlen(ret_str[3]) + 1, o);
     strcpy(ret_ptr, ret_str[3]);
     if (ret_str[3]) {
       free(ret_str[3]);
     }
     free(ret_str[0]);
     free(ret_str[1]);
     free(ret_str[2]);
   }
   free(duplex[0]);
   free(duplex[1]);
   free(duplex[2]);
   free(duplex[3]);

   return ret_ptr;
}
/* prints ascii output of hairpin structure */
std::string
drawHairpin(std::vector<int>& bp, double mh, double ms,
            const CProgParam_ThAl::mode mode, double temp, thal_results *o)
{
   int  ret_space = 0;
   std::string ret_str;
   int ret_last_l, ret_first_r, ret_center, ret_left_end, ret_right_start, ret_left_len, ret_right_len;
   int ret_add_sp_l, ret_add_sp_r;
   char ret_center_char;
   /* Plain text */
   int i, N;
   N = 0;
   double mg, t;
   if (!isFinite(ms) || !isFinite(mh)) {
      if((mode != CProgParam_ThAl::mode::FAST) && (mode != CProgParam_ThAl::mode::DEBUG_F))
      {
        if (mode != CProgParam_ThAl::mode::STRUCT)
        {
          printf("0\tdS = %g\tdH = %g\tinf\tinf\n", (double) ms,(double) mh);
#ifdef DEBUG
          fputs("No temperature could be calculated\n",stderr);
#endif
        }
      } else
      {
         o->temp = 0.0;       /* lets use generalization here */
         strcpy(o->msg, "No predicted sec struc for given seq\n");
      }
   } else
   {
      if((mode != CProgParam_ThAl::mode::FAST) && (mode != CProgParam_ThAl::mode::DEBUG_F))
      {
         for (i = 1; i < len1; ++i)
         {
            if(bp[i-1] > 0) N++;
         }
      } else {
         for (i = 1; i < len1; ++i)     // both branch equal ???
         {
            if(bp[i-1] > 0) N++;
         }
      }
      t = (mh / (ms + (((N/2)-1) * saltCorrection))) - ABSOLUTE_ZERO;

      if((mode != CProgParam_ThAl::mode::FAST) &&
         (mode != CProgParam_ThAl::mode::DEBUG_F))
      {
         mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));
         ms = ms + (((N/2)-1) * saltCorrection);
         o->temp = (double) t;

         if (mode != CProgParam_ThAl::mode::STRUCT)
         {
           printf("Calculated thermodynamical parameters for dimer:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
                  len1, (double) ms, (double) mh, (double) mg, (double) t);
         } else
         {
           sprintf(ret_para, "t: %.1f  dG: %.0f  dH: %.0f  dS: %.0f\\n",
                   (double) t, (double) mg, (double) mh, (double) ms);
         }
      }
      else
      {
         o->temp = (double) t;
         return {};
      }
   }
   /* plain-text output */
   std::string asciiRow (len1, '0');
   for(i = 1; i < len1+1; ++i) {
      if(bp[i-1] == 0) {
         asciiRow[(i-1)] = '-';
      } else {
         if(bp[i-1] > (i-1)) {
            asciiRow[(bp[i-1]-1)]='\\';
         } else  {
            asciiRow[(bp[i-1]-1)]='/';
         }
      }
   }
   if ((mode == CProgParam_ThAl::mode::GENERAL) || (mode == CProgParam_ThAl::mode::DEBUG))
   {
     printf("SEQ\t");
     for(i = 0; i < len1; ++i) printf("%c",asciiRow[i]);   // ??
     printf("\nSTR\t%s\n", oligo1);
   }
   if (mode == CProgParam_ThAl::mode::STRUCT) {
     ret_str = {};

     //save_append_string(&ret_str, &ret_space, o, ret_para);

     ret_last_l = -1;
     ret_first_r = -1;
     ret_center_char = '|';
     for(i = 0; i < len1; ++i)
     {
       if ( asciiRow[i] == '/')          ret_last_l  = i;
       if ((ret_first_r == -1 ) &&
           (asciiRow[i] == '\\'))        ret_first_r = i;
     }
     ret_center = ret_first_r - ret_last_l;

     if (ret_center % 2 == 0)        /* ret_center is odd */
     {
       ret_left_end    = ret_last_l + (ret_first_r - ret_last_l) / 2 - 1;
       ret_center_char = (char) oligo1[ret_left_end + 1]; 
       ret_right_start = ret_left_end + 2;
     }
     else
     {                                 /* ret_center is even */
       ret_left_end    = ret_last_l + (ret_first_r - ret_last_l - 1) / 2;
       ret_right_start = ret_left_end + 1;
     }

     ret_left_len  = ret_left_end + 1;
     ret_right_len = len1 - ret_right_start;
     ret_add_sp_l  = 0;
     ret_add_sp_r  = 0;

     if (ret_left_len > ret_right_len)       ret_add_sp_r = ret_left_len  - ret_right_len + 1;
     if (ret_right_len > ret_left_len)       ret_add_sp_l = ret_right_len - ret_left_len;

     for (i = 0 ; i < ret_add_sp_l ; i++)  ret_str += ' ' ; // save_append_char(&ret_str, &ret_space, o, ' ');
     ret_str += "5' " ;                                     // save_append_string(&ret_str, &ret_space, o, "5' ");
     for (i = 0 ; i < ret_left_len ; i++)  ret_str += (char) oligo1[i] ;
     ret_str += "U+2510\\n   " ;
     for (i = 0 ; i < ret_add_sp_l ; i++)  ret_str += ' ' ;
     for (i = 0 ; i < ret_left_len ; i++)
     {
       if (asciiRow[i] == '/')             ret_str += '|' ;
       else                                ret_str += ' ' ;
     }
     if (ret_center_char == '|' )          ret_str += "U+2502" ;
     else                                  ret_str += ret_center_char ;
     ret_str += "\\n" ;
     for (i = 0 ; i < ret_add_sp_r - 1 ; i++)       ret_str += ' ' ;
     ret_str += "3' " ;
     for (i = len1 ; i > ret_right_start - 1; i--)  ret_str += (char) oligo1[i] ;
     ret_str += "U+2518\\n" ;
/*
     save_append_string(&ret_str, &ret_space, o, "SEQ ");
     for(i = 0; i < len1; ++i) {
       save_append_char(&ret_str, &ret_space, o, asciiRow[i]);
     }
     save_append_string(&ret_str, &ret_space, o, "\\nSTR ");
     save_append_string(&ret_str, &ret_space, o, (const char*) oligo1);
     save_append_string(&ret_str, &ret_space, o, "\\n");
*/
   }
   return ret_str;
}

static void
save_append_string(char** ret, int *space, thal_results *o, const char *str) {
  int xlen, slen;
  if (str == NULL) {
    return;
  }
  if (*ret == NULL) {
    *ret = (char *) safe_malloc(sizeof(char)*500, o);
    *ret[0] = '\0';
    *space = 500;
  }
  xlen = strlen(*ret);
  slen = strlen(str);
  if (xlen + slen + 1 > *space) {
    *space += 4 * (slen + 1);
    *ret = (char *) safe_realloc(*ret, *space, o);
  }
  strcpy(*ret + xlen, str);
  return;
}

static void
save_append_char(char** ret, int *space, thal_results *o, const char str) {
  char fix[3];
  fix[0] = str;
  fix[1] = '\0';
  save_append_string(ret, space, o, fix);
}


static int 
equal(double a, double b)
{
#ifdef INTEGER
   return a == b;
#endif

   if (!isfinite(a) || !isfinite(b))
     return 0;
   return fabs(a - b) < 1e-5;

   if (a == 0 && b == 0)
     return 1;
}

static void 
strcatc(char* str, char c)
{
   str[strlen(str) + 1] = 0;
   str[strlen(str)] = c;
}


//     *******************    thal function  *****************************
/* central method: execute all sub-methods for calculating secondary
   structure for dimer or for monomer */
void thal( const seq& oligo_f,
           const seq& oligo_r,
           CProgParam_ThAl *a,
           const CProgParam_ThAl::mode modem,
           thal_results *o)
{
   if (!o) return; /* Leave it to the caller to crash */
   double* SH;
   int *bp;
   double mh, ms;
   double G1, bestG;

   enthalpyDPT = entropyDPT = NULL;
   numSeq1     = numSeq2    = {};
   oligo1      = oligo2     = {};

   o->msg = "";                          // ???
   o->temp = THAL_ERROR_SCORE;
   errno = 0;

    int len_f = oligo_f.length() ;
    int len_r = oligo_r.length() ;

   /* The following error messages will be seen by end users and will
      not be easy to understand. */
    if ((len_f > THAL_MAX_ALIGN) && (len_r > THAL_MAX_ALIGN))
      throw std::length_error ("Both sequences longer than " + std::itos( THAL_MAX_ALIGN) +
                               " for thermodynamic alignment");

    if (len_f > THAL_MAX_SEQ )
        throw std::length_error ("Target sequence (1) length > maximum allowed (" + std::itos( MAX_LEN)  ") "
                                 "in thermodynamic alignment"  );

    if (len_r > THAL_MAX_SEQ )
        throw std::length_error ("Target sequence (2) length > maximum allowed (" + std::itos( MAX_LEN)  ") "
                                 "in thermodynamic alignment"  );

    if (!a) throw std::runtime_error("Null 'in' argument pointer passed to thermodynamic alignment");

    if (  a->type != CProgParam_ThAl::type::Any      // 1
       && a->type != CProgParam_ThAl::type::end1     // 2
       && a->type != CProgParam_ThAl::type::end2     // 3
       && a->type != CProgParam_ThAl::type::Hairpin) // 4
          throw std::runtime_error("Illegal task type passed to thermodynamic alignment");

   o->align_end_1 = -1;
   o->align_end_2 = -1;

   if (!oligo_f ) throw std::invalid_argument("Empty first sequence passed to thermodynamic alignment");
   if (!oligo_r ) throw std::invalid_argument("Empty second sequence passed to thermodynamic alignment");

   oligo1 = oligo_f;
   oligo2 = oligo_r;
   if(a->type!=CProgParam_ThAl::type::end2) oligo1.swap(oligo2);     // 3
   len1 = oligo1.length();
   len2 = oligo2.length();
   for(int i = 0; i < len1; i++) oligo1[i] = toupper(oligo1[i]);
   for(int i = 0; i < len2; i++) oligo2[i] = toupper(oligo2[i]);

    /*** INIT values for unimolecular and bimolecular structures ***/

   if (a->type==CProgParam_ThAl::type::Hairpin)                         /* unimolecular folding     4 */
   {
      len3 = len2 -1;
      dplx_init_H = 0.0;                /* initiation enthalpy; for duplex 200, for unimolecular structure 0 */
      dplx_init_S = -0.00000000001;     /* initiation entropy; for duplex -5.7, for unimoleculat structure 0 */
      RC=0;
      send5.resize(len1);   hend5.resize(len1);
   }
   else //if(a->type!=4)                           /* all the others ----> hybridization of two oligos */
   {                                               /* using the second reversed   */
      reverse(oligo2);                  /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
      dplx_init_H = 200;
      dplx_init_S = -5.7;
      if(symmetry_thermo(oligo1) && symmetry_thermo(oligo2))
      {
         RC = R  * log(a->dna_conc/1000000000.0);      // both symmetric
      }
      else
      {
         RC = R  * log(a->dna_conc/4000000000.0);
      }
   }
   /*** Calc part of the salt correction ***/
   saltCorrection=saltCorrectS(a->mv,a->dv,a->dntp); /* salt correction for entropy, must be multiplied with N, which is
                                              the total number of phosphates in the duplex divided by 2; 8bp dplx N=7 */

   numSeq1 = '\4' + nt2code(oligo1) + '\4'    /* mark first and last as N-s */
   numSeq2 = '\4' + nt2code(oligo2) + '\4'    /* mark first and last as N-s */

   if (a->type == CProgParam_ThAl::type::Hairpin)   /* 4 calculate structure of monomer */
   {
      enthalpyDPT.clear(); enthalpyDPT.resize( len1 * len2 );   // safe_realloc(ptr, m * n * sizeof(double)
      entropyDPT.clear() ; entropyDPT.resize ( len1 * len2 );
      initMatrix2();                   // initiates thermodynamic parameter tables of entropy and enthalpy for monomer
      fillMatrix2(a->maxLoop);         // calc-s thermod values into dynamic progr table (monomer). MAJOR CALCULATIONS --> maxTM2(), CBI(), calc_bulge_internal2 !!
      calc_terminal_bp(a->temp);
      mh = HEND5(len1);
      ms = SEND5(len1);
      o->align_end_1 = (int) mh;
      o->align_end_2 = (int) ms;
      auto bp = std::vector<int> (len1, 0);
      if(isFinite(mh))                          //  RESULTS HERE   !!!!!!!!!!
      {
         tracebacku(bp, a->maxLoop, o);    // traceback for hairpins in unimolecular structure. calls  /* traceback for unimolecular structure */
         o->sec_struct=drawHairpin(bp, mh, ms, mode,a->temp, o); /* if mode=THL_FAST or THL_DEBUG_F then return after printing basic therm data */
          /* prints ascii output of hairpin structure */
      }
      else if((mode != CProgParam_ThAl::mode::FAST   ) &&
              (mode != CProgParam_ThAl::mode::DEBUG_F) &&
              (mode != CProgParam_ThAl::mode::STRUCT)    )
      {
         fputs("No secondary structure could be calculated\n",stderr);
      }

      if(o->temp==-_INFINITY && (!strcmp(o->msg, ""))) o->temp=0.0;

      enthalpyDPT={}; // .clear(); // shrink_to_fit   ??   free all memmory ?!
      entropyDPT ={}; // .clear(); // shrink_to_fit   ??
      numSeq1    ={}; // .clear(); // shrink_to_fit   ??
      numSeq2    ={}; // .clear(); // shrink_to_fit   ??
      send5      ={}; // .clear(); // shrink_to_fit   ??
      hend5      ={}; // .clear(); // shrink_to_fit   ??
      oligo1     ={}; // .clear(); // shrink_to_fit   ??
      oligo2     ={}; // .clear(); // shrink_to_fit   ??

      return;

   }
   else          // if(a->type!=4) {                             /* Hybridization of two moleculs */
   {
      len3 = len2;                 // safe_realloc(ptr, m * n * sizeof(double)
      enthalpyDPT.clear() ; enthalpyDPT.resize( len1 * len2 );  /* dyn. programming table for dS and dH */
      entropyDPT .clear() ; entropyDPT .resize( len1 * len2 );  /* enthalpyDPT is 3D array represented as 1D array */

      initMatrix();               ///< initiates thermodynamic parameter tables of entropy and enthalpy for dimer
      fillMatrix(a->maxLoop );    ///< initiates thermodynamic parameter tables of entropy and enthalpy for monomer
      double SH[2];
      /* calculate terminal basepairs */
      bestI = bestJ = 0;
      G1 = bestG = _INFINITY;
      if(a->type==1)
         for (i = 1; i <= len1; i++) {
            for (j = 1; j <= len2; j++) {
               RSH(i, j, SH);
               SH[0] = SH[0]+SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0 */
               SH[1] = SH[1]+SMALL_NON_ZERO;
               G1 = (EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) - TEMP_KELVIN*(EntropyDPT(i, j) + SH[0] + dplx_init_S);
               if(G1<bestG){
                  bestG = G1;
                  bestI = i;
                  bestJ = j;
               }
            }
         }
      std::vector<int> ps1(len1, 0), ps2(len2, 0);
      if(a->type == CProgParam_ThAl::type::end1 || a->type == CProgParam_ThAl::type::end2)        {
         /* THAL_END1 */
         bestI = bestJ = 0;
         bestI = len1;
         i = len1;
         G1 = bestG = _INFINITY;
         for (j = 1; j <= len2; ++j) {
            RSH(i, j, SH);
            SH[0] = SH[0]+SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0,
                                             that compiler could understand that SH is changed in this cycle */
            SH[1] = SH[1]+SMALL_NON_ZERO;
            G1 = (EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) - TEMP_KELVIN*(EntropyDPT(i, j) + SH[0] + dplx_init_S);
            if(G1<bestG){
               bestG = G1;
               bestJ = j;
            }
         }
      }
      if (!isFinite(bestG)) bestI = bestJ = 1;
      double dH, dS;
      RSH(bestI, bestJ, SH);
      dH =  EnthalpyDPT(bestI, bestJ)+ SH[1] + dplx_init_H;
      dS = (EntropyDPT (bestI, bestJ)+ SH[0] + dplx_init_S);

      for (i = 0; i < len1; ++i)      ps1[i] = 0;           /* tracebacking */
      for (j = 0; j < len2; ++j)      ps2[j] = 0;

      if(isFinite(EnthalpyDPT(bestI, bestJ)))                         //  RESULTS HERE   !!!!!!!!!!
      {
         traceback(bestI, bestJ, RC, ps1, ps2, a->maxLoop, o);   /* traceback for dimers */
         o->sec_struct=drawDimer(ps1, ps2, SHleft, dH, dS, mode, a->temp, o);
         o->align_end_1=bestI;
         o->align_end_2=bestJ;
      } else  {
         o->temp = 0.0;
         /* fputs("No secondary structure could be calculated\n",stderr); */
      }
      free(ps1);
      free(ps2);
      free(SH);
      free(enthalpyDPT);
      free(entropyDPT);
      numSeq1.clear();
      numSeq2.clear();
      oligo1.clear();
      return;
   }
   return;
}
/*** END thal() ***/


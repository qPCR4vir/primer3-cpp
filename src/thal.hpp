/* Copyright (c) 1996 - 2018
 Whitehead Institute for Biomedical Research, Steve Rozen, Andreas Untergasser
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
 All rights reserved.

     This file is part the primer3 software suite.

     This software suite is free software;
     you can redistribute is and/or modify it under the terms
     of the GNU General Public License as published by the Free
     Software Foundation; either version 2 of the License, or (at
     your option) any later version.

     This software is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this file (file gpl-2.0.txt in the source
     distribution); if not, write to the Free Software
     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, ST
 RICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#ifndef _THAL_H
#define _THAL_H
#include <string>
#include <istream>
#include <memory>

#include "filesystem.hpp"

#include <float.h> /* ! mul ei ole float.h-d includes DBL_MAX */
#include <math.h>
#include <limits.h>

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

#ifndef THAL_ERROR_SCORE
# define THAL_ERROR_SCORE -_INFINITY
#endif

/* The maximum length of _one_ of the two sequences being aligned in a
   thermodynamic alignment. In other words, the length of one sequence
   must be <= THAL_MAX_ALIGN, but the other sequence can be longer.
   The rationale behind this value (60) is that this is the maxium
   reasonable length for nearest neighbor models. It is the maxium
   length at which we can restrict our model to only two states of
   melting: fully intact duplex or completely dissociated single
   strands. */
#ifndef THAL_MAX_ALIGN
#define THAL_MAX_ALIGN 60
#endif

/* The maxium length of the other sequence in a thermodynamic
   alignment. This value can be increased, though alignments against
   very long sequences will be quite slow. As of 2012-05-18, we only
   potentially see sequences longer this when checking for mispriming
   in the template ('max_template_mispriming') in libprimer3.c, which
   is really designed to find sites of ectopic primer very close (a
   few kilobases) from the location of the cadidate primer. */
#ifndef THAL_MAX_SEQ
#define THAL_MAX_SEQ   10000
#endif


extern const double ABSOLUTE_ZERO;
extern const int MIN_LOOP;
constexpr int    MAX_LOOP = 30; ///< the maximum size of loop that can be calculated; for larger loops formula must be implemented
constexpr double TEMP_KELVIN = 310.15;

/// handle the NN parameters for the thermodynamic alignment
class thal_parameters
 {
public:
    thal_parameters();                                                 ///< calls set_defaults( );
    explicit thal_parameters(const std::filesystem::path& dirname );  ///< calls load( );
    ~thal_parameters();

    /// set hard coded defaults
    int set_defaults( );

    /// load and parse parameters from a directory
    int load(const std::filesystem::path& dirname );

    class impl;
    std::unique_ptr<impl> m_thal_p;
} ;

/// for compatibility
void set_thal_default_args(thal_parameters *a)
{
    new (a) thal_parameters;
}
//void set_thal_oligo_default_args(thal_args *a);


/// Structure for passing arguments to THermodynamic ALignment calculation
class thal_args
{

 public:
    const thal_parameters& tp;
    enum class type {Any =1, end1, end2, Hairpin };

   type         type = type::Any;  ///< 1 = Any, (by default)
   int          maxLoop{MAX_LOOP}; ///< maximum size of loop to consider; longer than 30 bp are not allowed
   double       mv   = 50;         ///< mM, concentration of monovalent cations
   double       dv   = 0.0;        ///< mM, cconcentration of divalent cations
   double       dntp = 0.8;        ///< mM, cconcentration of dNTP-s
   double       dna_conc = 50;     ///< nM, cconcentration of oligonucleotides
   double       temp = TEMP_KELVIN;///< Kelvin, ctemperature from which hairpin structures will be calculated
   int          dimer = 1;         ///< if non zero, dimer structure is calculated

    enum class mode {
        FAST    = 0, //< = 0 - score only with optimized functions (fast)
        GENERAL = 1, //< = 1 - use general function without debug (slow)
        DEBUG_F = 2, //< = 2 - debug mode with fast, print alignments on STDERR
        DEBUG_  = 3, //< = 3 - debug mode print alignments on STDERR
        STRUCT  = 4  //< = 4 - calculate secondary structures as string
    };

    explicit thal_args (const thal_parameters& tp) : tp(tp){};

    void set_defaults      ( )
    {
        this->type     = type::Any;     /* thal_alignment_type THAL_ANY */
        this->maxLoop  = MAX_LOOP;
        this->mv       = 50;            /* mM */
        this->dv       = 0.0;           /* mM */
        this->dntp     = 0.8;           /* mM */
        this->dna_conc = 50;            /* nM */
        this->temp     = TEMP_KELVIN;   /* Kelvin */
        this->dimer    = 1;             /* by default dimer structure is calculated */
    }

    void set_oligo_defaults( )
    {
        this->type     = type::Any;     /* thal_alignment_type THAL_ANY */
        this->maxLoop  = MAX_LOOP;
        this->mv       = 50;            /* mM */
        this->dv       = 0.0;           /* mM */
        this->dntp     = 0.0;           /* mM the only difference !!!! */
        this->dna_conc = 50;            /* nM */
        this->temp     = TEMP_KELVIN;   /* Kelvin */
        this->dimer    = 1;             /* by default dimer structure is calculated */
    }


    ~thal_args() = default;
    // int  thal_free_parameters(thal_parameters *a);
} ;

/* Structure for receiving results from the thermodynamic alignment calculation */
class thal_results
{
public:
    std::string  msg;
    double       temp;
    int          align_end_1;
    int          align_end_2;
    std::string  sec_struct;
} ;

using seq = std::basic_string<unsigned char> ;

///  Central method for finding the best alignment.  On error, o->temp
///    is set to THAL_ERROR_SCORE and a message is put in o->msg.  The
///    error might be caused by ENOMEM. To determine this it is necessary
///    to check errno.
void thal( const seq& oligo_f,
           const seq& oligo_r,
           const thal_args &a,
           thal_args::mode mode,
           thal_results &o);

#endif

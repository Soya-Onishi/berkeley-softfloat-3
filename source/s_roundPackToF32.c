
/*============================================================================

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3e, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2017 The Regents of the University of
California.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions, and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions, and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the University nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=============================================================================*/

#include <stdbool.h>
#include <stdint.h>
#include "platform.h"
#include "internals.h"
#include "softfloat.h"

float32_t
 softfloat_roundPackToF32( bool sign, int_fast16_t exp, uint_fast32_t sig )
{
    uint_fast8_t roundingMode;
    bool roundNearEven;
    uint_fast8_t roundIncrement, roundBits;
    bool isTiny;
    uint_fast32_t uiZ;
    union ui32_f32 uZ;

    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
    roundingMode = softfloat_roundingMode;
    roundNearEven = (roundingMode == softfloat_round_near_even);
    roundIncrement = 0x40;

    // 最近接丸めでは無い時
    if ( ! roundNearEven && (roundingMode != softfloat_round_near_maxMag) ) {
        /*
          | mode | sign | roundIncrement | 
          |  min |    + |           0x00 |
          |  min |    - |           0x7F |
          |  max |    + |           0x7F |
          |  max |    - |           0x00 |
          | zero |  +/- |           0x00 |
        */        
        uint_fast8_t mode = sign ? softfloat_round_min : softfloat_round_max;
        roundIncrement = mode == roundingMode ? 0x7F : 0;
    }

    // softfloat_addMagsF32 ではこちらに渡される際に sig の値は0b01xx_xxxx_....となっている
    // つまり仮数部はいかのようになる
    // 0b01 | 仮数部(23bit) | xxx xxxx 
    // よって下位7ビットが丸めビットとなる．    
    roundBits = sig & 0x7F;
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/

   // 0xFD以上になるのは exp が 0xFD, 0xFE, および負の値のときのみ
   // exp は 0xFF にはならない．
   // なぜなら，この関数を呼ぶ前にそのような値はif文によってフィルタリングされるから
   // (加算の話，乗算とかだとexpが0xFFになったりするかも)
    if ( 0xFD <= (unsigned int) exp ) {
        if ( exp < 0 ) {
            /*----------------------------------------------------------------
            *----------------------------------------------------------------*/
            isTiny =
                (softfloat_detectTininess == softfloat_tininess_beforeRounding)
                    || (exp < -1) || (sig + roundIncrement < 0x80000000);
            sig = softfloat_shiftRightJam32( sig, -exp );
            exp = 0;
            roundBits = sig & 0x7F;
            if ( isTiny && roundBits ) {
                softfloat_raiseFlags( softfloat_flag_underflow );
            }
        } else if ( (0xFD < exp) || (0x80000000 <= sig + roundIncrement) ) {
            /*----------------------------------------------------------------
            *----------------------------------------------------------------*/
            // 0xFEでも強制的にオーバーフロー判定になるのは sig の値がすでに繰り上がった状態になっているから．
                        
            softfloat_raiseFlags(
                softfloat_flag_overflow | softfloat_flag_inexact );
            
            // オーバーフロー時に返す値は丸めモードや符号によって異なる
            // https://docs.oracle.com/cd/E19957-01/806-4847/ncg_handle.html
            // 
            // | roundIncrement |             uiZ            |
            // |              0 | | sign | 0xFE | 0x7FFFFF | |
            // |      otherwise | | sign | 0xFF | 0x000000 | |
            uiZ = packToF32UI( sign, 0xFF, 0 ) - ! roundIncrement;
            goto uiZ;
        }
    }
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
    
    sig = (sig + roundIncrement)>>7;

    // 丸めが発生したら不正確例外フラグをセットする
    if ( roundBits ) {
        softfloat_exceptionFlags |= softfloat_flag_inexact;
#ifdef SOFTFLOAT_ROUND_ODD
        if ( roundingMode == softfloat_round_odd ) {
            sig |= 1;
            goto packReturn;
        }
#endif
    }
    
    // 最近接丸めの際にLSBを0にするための処理
    sig &= ~(uint_fast32_t) (! (roundBits ^ 0x40) & roundNearEven);

    // 謎
    if ( ! sig ) exp = 0;
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
 packReturn:
    // sig のケチ表現が残ったままだけど，計算によってケチ表現の桁上がりも packToF32UI 上で吸収してくれる
    // 言いたいことが書けないけど，正しく計算できているということだけは確か．（多少のトリッキーさは感じるけど）
    uiZ = packToF32UI( sign, exp, sig );
 uiZ:
    uZ.ui = uiZ;
    return uZ.f;

}



/*============================================================================

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3e, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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
#include "specialize.h"

float32_t softfloat_addMagsF32( uint_fast32_t uiA, uint_fast32_t uiB )
{
    int_fast16_t expA;
    uint_fast32_t sigA;
    int_fast16_t expB;
    uint_fast32_t sigB;
    int_fast16_t expDiff;
    uint_fast32_t uiZ;
    bool signZ;
    int_fast16_t expZ;
    uint_fast32_t sigZ;
    union ui32_f32 uZ;

    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
    expA = expF32UI( uiA );
    sigA = fracF32UI( uiA );
    expB = expF32UI( uiB );
    sigB = fracF32UI( uiB );
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
    expDiff = expA - expB;

    // expA と expB が同じ時
    if ( ! expDiff ) {
        /*--------------------------------------------------------------------
        *--------------------------------------------------------------------*/

       // expA と expB が 0x00 のとき
       // sigA が0のときは => 0
       //       それ以外は => 非正規化数
        if ( ! expA ) {
            // 非正規化数は足すだけでいいのでこれでいい．
            // 仮に桁が上がったとしても
            // 指数部が桁上がりによって 0000_0001 になるので正しい．
            uiZ = uiA + sigB;
            goto uiZ;
        }

        // expA と expB が 0xFF のとき
        // sigA または sigB が0以外ならNaN
        // sigA と sigB が0なら無限大
        if ( expA == 0xFF ) {
            if ( sigA | sigB ) goto propagateNaN;
            uiZ = uiA;
            goto uiZ;
        }

        // 以降は expA と expB が0x00でも0xFFでもない場合
        signZ = signF32UI( uiA );
        expZ = expA;

        // ケチ表現によって消えている 1.xx の部分が足されたことで 2.xx になった
        // そのため，0x0100_0000 によって 25ビット目に1をセットしている（ケチ表現の1.xxは24ビット目に1をセットしていることに相当する）
        sigZ = 0x01000000 + sigA + sigB;

        // !(sigZ & 1) => 奇数じゃないか見てる？
        // exp < 0xFE  => なんで exp <= 0xFE じゃないの？
        if ( ! (sigZ & 1) && (expZ < 0xFE) ) {
            uiZ = packToF32UI( signZ, expZ, sigZ>>1 );
            goto uiZ;
        }
        
        sigZ <<= 6;
    } else {
        /*--------------------------------------------------------------------
        *--------------------------------------------------------------------*/
        signZ = signF32UI( uiA );

        // 丸めビット（ガードビットなど）の3ビット分シフトさせるのは分かるけど
        // なぜ6ビットもシフトさせるのか
        //
        // 6ビットもシフトさせると softfloat_shiftRightJam32 内で右シフトにより消えたビットを検査して
        // スティッキービットと思われる部分に1をセットする操作がうまく行かない場合があるのでないか
        // 例:
        //   右に4ビットシフト
        //   001 | 23bit | 000000 => 0000 | 001 | 19bit | xxxx00
        //   xxxx が 0001 であった場合，3ビットのみのシフトであればスティッキービットに1がセットされた
        //   しかし，6ビットもシフトさせたために，スティッキービットに1をセットされなくなった．        
        sigA <<= 6;
        sigB <<= 6;
        if ( expDiff < 0 ) {
            // sigA のほうが小さいので，sigA を右シフトする

            if ( expB == 0xFF ) {
                if ( sigB ) goto propagateNaN;
                uiZ = packToF32UI( signZ, 0xFF, 0 );
                goto uiZ;
            }
            expZ = expB;

            // 非正規化数の場合は仮数部を2倍
            // おそらく softfloat_roundPackToF32 で下位7bitが丸めビットとして扱われるので，それの対処
            // (2倍は左に1ビットシフトしたことと同義であることに注意)
            // 
            // それ以外のときは仮数部の24ビット目に1をセット（ケチ表現の解消）
            // 上で sigA は 6ビットシフトされていることに注意
            // ここでは0x2000_0000を足すのみで7bitの対処をしてないが，下の sigZ < 0x4000_0000 の部分で行われる（と，思う）．
            sigA += expA ? 0x20000000 : sigA;
            sigA = softfloat_shiftRightJam32( sigA, -expDiff );
        } else {
            if ( expA == 0xFF ) {
                if ( sigA ) goto propagateNaN;
                uiZ = uiA;
                goto uiZ;
            }
            expZ = expA;
            sigB += expB ? 0x20000000 : sigB;
            sigB = softfloat_shiftRightJam32( sigB, expDiff );
        }

        // 右シフトしなかった方の仮数部のケチ表現分も加算（0x2000_0000）
        sigZ = 0x20000000 + sigA + sigB;
        if ( sigZ < 0x40000000 ) {
            // 非正規化数の場合は 0x4000_0000 にはならないので
            // sigZ を左シフトして指数部の値を1減らす．
            --expZ;
            sigZ <<= 1;
        }
    }
    return softfloat_roundPackToF32( signZ, expZ, sigZ );
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
 propagateNaN:
    uiZ = softfloat_propagateNaNF32UI( uiA, uiB );
 uiZ:
    uZ.ui = uiZ;
    return uZ.f;

}


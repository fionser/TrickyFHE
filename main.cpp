#include <NTL/ZZX.h>
#include "HElib/FHE.h"

/// Symmetric encryption for m.
/// Basically using a RLWE instance (c0, c1) where c0 + c1 * s = q e for a short e
/// Then, the ciphertext of m is (c0 + m, c1) (over the ciphertext space Z[X]_Q / (X^n + 1) for a large Q)
std::vector<NTL::ZZX> encrypt(NTL::ZZX const& m,
                              DoubleCRT const& s,
                              FHEcontext const& context) {
    std::vector<CtxtPart> parts;
    {
        CtxtPart tmp(context, context.ctxtPrimes);
        parts.assign(2, tmp);
    }

    long ptxtSpace = context.alMod.getPPowR();
    RLWE(parts[0], parts[1], s, ptxtSpace);
    parts[0] += m;

    std::vector<NTL::ZZX> ret(parts.size(), NTL::to_ZZX(0));
    for (size_t i = 0; i < parts.size(); i++)
        // should use false here, because the error terms are taken from gaussian
        parts[i].toPoly(ret[i], context.ctxtPrimes, false);

    return ret;
}

// The ciphertext encrypts a polynomial M, to create a ciphertext that encrypts the specific coefficients M[loc].
NTL::ZZX extract(std::vector<NTL::ZZX> const& parts, 
                 long loc,
                 FHEcontext const& context) {
    if (parts.size() != 2)
        return NTL::to_ZZX(0);
    long phim = context.zMStar.getPhiM();
    assert(loc >= 0 && loc < phim);

    NTL::ZZX ret;
    ret.SetLength(phim + 1);
    NTL::SetCoeff(ret, 0, NTL::coeff(parts[1], loc));

    // loc = i + j \mod phim
    // j = loc - i \mod phim
    for (long i = 0; i < phim; i++) {
        long j = loc - i;
        long jj = j < 0 ? j + phim : j;
        auto tmp = NTL::coeff(parts[1], jj % phim);
        if (j < 0 || j >= phim)
            NTL::negate(tmp, tmp);
        NTL::SetCoeff(ret, i, tmp);
    }
    NTL::SetCoeff(ret, phim, NTL::coeff(parts[0], loc));
    return ret;
}

// Decrypt the extracted ciphertext.
long decrypt(NTL::ZZX const& ctx,
             NTL::ZZX const& s,
             FHEcontext const& context) {
    long phim = context.zMStar.getPhiM();
    long ptxtSpace = context.alMod.getPPowR();
    auto Q = context.productOfPrimes(context.ctxtPrimes);
    auto halfQ = Q >> 1;
    
    NTL::ZZ inner_product(NTL::coeff(ctx, phim));
    for (long i = 0; i < phim; i++) {
        inner_product += ctx[i] * s[i];
    }
    inner_product %= Q; // The operation is done in [-Q/2, Q/2)
    if (inner_product >= halfQ)
        inner_product -= Q;
    return NTL::rem(inner_product, ptxtSpace);
}

int main() {
    long m = 64;
    long phim = m >> 1;
    FHEcontext context(m, 101, 1);
    buildModChain(context, 3);

    DoubleCRT s(context);
    s.sampleHWt(16);
    NTL::ZZX s_poly;
    s.toPoly(s_poly, false);

    FHESecKey sk(context);
    sk.ImportSecKey(s, 64);

    NTL::ZZX mess, mess2;
    mess.SetLength(phim);
    for (long i = 0; i < phim; i++)
        mess[i] = NTL::RandomBnd(101);

    auto parts = encrypt(mess, s, context);
    for (long i = 0; i < phim; i++) {
        auto ctx = extract(parts, i, context);
        if (decrypt(ctx, s_poly, context) != mess[i]) 
            std::cout << "fail at " << i << " loc\n";
    }
    return 0;
}

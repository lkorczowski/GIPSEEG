function SYNC=SYNC_Riem(P1,P2)
        Cref=blkdiag(cov(P1'),cov(P2'));
        C=cov([P1;P2]');
        SYNC=distance(Cref,C,'riemann');
end
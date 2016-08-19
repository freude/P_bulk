function M=make_pot_mat(j1,j2,jj1,jj2,jjj1,jjj2,amp,MEs)

    if (j1~=j2)
        M=MEs(j1,j2)*amp(jjj1,jj1)*amp(jjj2,jj2);
    end;

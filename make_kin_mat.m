function M=make_kin_mat(j1,j2,jj1,jj2,jjj1,jjj2,Energy,A,B,S)

    if (j1==j2)
        if ((j2==jjj1)&&(jj1==jj2)&&(jjj1==jjj2))
            M=Energy(jj1);
        end;

        if ((j2==jjj1)&&(jj1~=jj2)&&(jjj1==jjj2))
            M=0;
        end;

        if ((j2~=jjj1)&&(j2~=jjj2))
            if (jjj1==jjj2)
                M=A(jj1,jj2);
            else
                M=B(jj1,jj2);
            end;
        end;

        if ((j2==jjj1)&&(j2~=jjj2))
            M=Energy(jj1)*S(jjj1,jjj2,jj1,jj2);
        end;

        if ((j2==jjj2)&&(j2~=jjj1))
            M=Energy(jj2)*S(jjj1,jjj2,jj1,jj2);
        end;
    end;

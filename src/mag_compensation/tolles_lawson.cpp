/***********************************************
#
# Author: Sun Qinxuan
#
# Email: sunqinxuan@outlook.com
#
# Last modified:	2023-08-15 14:04
#
# Filename:		tolles_lawson.cpp
#
# Description:
#
************************************************/

#include "mag_compensation/tolles_lawson.hpp"

namespace magnav {

bool TollesLawson::createMatrixA(const std::vector<double> &Bx,const std::vector<double> &By,const std::vector<double> &Bz,
																 std::vector<double> &Bt,
																 const std::unordered_set<TLterm> &terms,
																 const double Bt_scale,
																 const bool return_B)
{
	if(!(Bx.size()==By.size() && Bx.size()==Bz.size()))
	{
		ERROR("[TollesLawson][createMatrixA] Bx By Bz sizes not equal!");
		return false;
	}
	if(terms.empty())
	{
		ERROR("[TollesLawson][createMatrixA] terms empty!");
		return false;
	}

	if(Bt.empty())
	{
		Bt.resize(Bx.size());
		for(size_t i=0;i<Bx.size();i++)
		{
			Bt[i]=sqrt(Bx[i]*Bx[i]+By[i]*By[i]+Bz[i]*Bz[i]);
		}
	}
	else
	{
		if(Bx.size()!=Bt.size())
		{
			ERROR("[TollesLawson][createMatrixA] Bx Bt sizes not equal!");
			return false;
		}
	}

	int N=Bx.size();
	std::vector<double> Bx_hat(N), By_hat(N), Bz_hat(N);
	for(int i=0;i<N;i++)
	{
		Bx_hat[i]=Bx[i]/Bt[i];
		By_hat[i]=By[i]/Bt[i];
		Bz_hat[i]=Bz[i]/Bt[i];
	}

	std::vector<double> Bx_dot,By_dot,Bz_dot;
	fdm(Bx_dot,Bx);
	fdm(By_dot,By);
	fdm(Bz_dot,Bz);

	std::vector<double>  Bx_hat_Bx(N), Bx_hat_By(N), Bx_hat_Bz(N), By_hat_By(N), By_hat_Bz(N), Bz_hat_Bz(N),
    Bx_hat_Bx_dot(N), Bx_hat_By_dot(N), Bx_hat_Bz_dot(N), By_hat_Bx_dot(N), By_hat_By_dot(N), By_hat_Bz_dot(N), Bz_hat_Bx_dot(N), Bz_hat_By_dot(N), Bz_hat_Bz_dot(N);
	for(int i=0;i<N;i++)
	{
    Bx_hat_Bx[i] = Bx_hat[i] * Bx[i] / Bt_scale;
    Bx_hat_By[i] = Bx_hat[i] * By[i] / Bt_scale;
    Bx_hat_Bz[i] = Bx_hat[i] * Bz[i] / Bt_scale;
    By_hat_By[i] = By_hat[i] * By[i] / Bt_scale;
    By_hat_Bz[i] = By_hat[i] * Bz[i] / Bt_scale;
    Bz_hat_Bz[i] = Bz_hat[i] * Bz[i] / Bt_scale;

    Bx_hat_Bx_dot[i] = Bx_hat[i] * Bx_dot[i] / Bt_scale;
    Bx_hat_By_dot[i] = Bx_hat[i] * By_dot[i] / Bt_scale;
    Bx_hat_Bz_dot[i] = Bx_hat[i] * Bz_dot[i] / Bt_scale;
    By_hat_Bx_dot[i] = By_hat[i] * Bx_dot[i] / Bt_scale;
    By_hat_By_dot[i] = By_hat[i] * By_dot[i] / Bt_scale;
    By_hat_Bz_dot[i] = By_hat[i] * Bz_dot[i] / Bt_scale;
    Bz_hat_Bx_dot[i] = Bz_hat[i] * Bx_dot[i] / Bt_scale;
    Bz_hat_By_dot[i] = Bz_hat[i] * By_dot[i] / Bt_scale;
    Bz_hat_Bz_dot[i] = Bz_hat[i] * Bz_dot[i] / Bt_scale;
	}

	/*
void Eigen::PlainObjectBase< Derived >::resize	(	Index 	rows, Index 	cols )		inline

Resizes *this to a rows x cols matrix.

This method is intended for dynamic-size matrices, 
although it is legal to call it on any matrix as long as fixed dimensions are left unchanged. 
If you only want to change the number of rows and/or of columns, 
you can use resize(NoChange_t, Index), resize(Index, NoChange_t).

!!!
--------------------------------------------------------------------------------------------------
If the current number of coefficients of *this exactly matches the product rows * cols, 
then no memory allocation is performed and the current values are left unchanged. 
In all other cases, including shrinking, the data is reallocated and all previous values are lost.
--------------------------------------------------------------------------------------------------
*/
	Eigen::MatrixXd A;

	if(terms.find(PERMANENT)!=terms.end())
	{
		A.resize(N,3);
	}

	return true;
	/*
function create_TL_A(Bx, By, Bz;
                     Bt       = sqrt.(Bx.^2+By.^2+Bz.^2),
                     terms    = [:permanent,:induced,:eddy],
                     Bt_scale = 50000,
                     return_B = false)

    Bx_hat = Bx ./ Bt
    By_hat = By ./ Bt
    Bz_hat = Bz ./ Bt

    Bx_dot = fdm(Bx)
    By_dot = fdm(By)
    Bz_dot = fdm(Bz)

    Bx_hat_Bx = Bx_hat .* Bx ./ Bt_scale
    Bx_hat_By = Bx_hat .* By ./ Bt_scale
    Bx_hat_Bz = Bx_hat .* Bz ./ Bt_scale
    By_hat_By = By_hat .* By ./ Bt_scale
    By_hat_Bz = By_hat .* Bz ./ Bt_scale
    Bz_hat_Bz = Bz_hat .* Bz ./ Bt_scale

    Bx_hat_Bx_dot = Bx_hat .* Bx_dot ./ Bt_scale
    Bx_hat_By_dot = Bx_hat .* By_dot ./ Bt_scale
    Bx_hat_Bz_dot = Bx_hat .* Bz_dot ./ Bt_scale
    By_hat_Bx_dot = By_hat .* Bx_dot ./ Bt_scale
    By_hat_By_dot = By_hat .* By_dot ./ Bt_scale
    By_hat_Bz_dot = By_hat .* Bz_dot ./ Bt_scale
    Bz_hat_Bx_dot = Bz_hat .* Bx_dot ./ Bt_scale
    Bz_hat_By_dot = Bz_hat .* By_dot ./ Bt_scale
    Bz_hat_Bz_dot = Bz_hat .* Bz_dot ./ Bt_scale

    # # original (slightly incorrect) eddy current terms
    # Bx_hat_Bx_dot = Bx_hat .* fdm(Bx_hat) .* Bt ./ Bt_scale
    # Bx_hat_By_dot = Bx_hat .* fdm(By_hat) .* Bt ./ Bt_scale
    # Bx_hat_Bz_dot = Bx_hat .* fdm(Bz_hat) .* Bt ./ Bt_scale
    # By_hat_Bx_dot = By_hat .* fdm(Bx_hat) .* Bt ./ Bt_scale
    # By_hat_By_dot = By_hat .* fdm(By_hat) .* Bt ./ Bt_scale
    # By_hat_Bz_dot = By_hat .* fdm(Bz_hat) .* Bt ./ Bt_scale
    # Bz_hat_Bx_dot = Bz_hat .* fdm(Bx_hat) .* Bt ./ Bt_scale
    # Bz_hat_By_dot = Bz_hat .* fdm(By_hat) .* Bt ./ Bt_scale
    # Bz_hat_Bz_dot = Bz_hat .* fdm(Bz_hat) .* Bt ./ Bt_scale

    A = Array{eltype(Bt)}(undef,size(Bt,1),0)

    # add (3) permanent field terms - all
    if any([:permanent,:p,:permanent3,:p3] .∈ (terms,))
    	A = [A Bx_hat By_hat Bz_hat]
    end

    # add (6) induced field terms - all
    if any([:induced,:i,:induced6,:i6] .∈ (terms,))
        A = [A Bx_hat_Bx Bx_hat_By Bx_hat_Bz By_hat_By By_hat_Bz Bz_hat_Bz]
    end

    # add (5) induced field terms - all except Bz_hat_Bz
    if any([:induced5,:i5] .∈ (terms,))
        A = [A Bx_hat_Bx Bx_hat_By Bx_hat_Bz By_hat_By By_hat_Bz]
    end

    # add (3) induced field terms - Bx_hat_Bx, By_hat_By, Bz_hat_Bz
    if any([:induced3,:i3] .∈ (terms,))
        A = [A Bx_hat_Bx By_hat_By Bz_hat_Bz]
    end

    # add (9) eddy current terms - all
    if any([:eddy,:e,:eddy9,:e9] .∈ (terms,))
        A = [A Bx_hat_Bx_dot Bx_hat_By_dot Bx_hat_Bz_dot]
        A = [A By_hat_Bx_dot By_hat_By_dot By_hat_Bz_dot]
        A = [A Bz_hat_Bx_dot Bz_hat_By_dot Bz_hat_Bz_dot]
    end

    # add (8) eddy current terms - all except Bz_hat_Bz_dot
    if any([:eddy8,:e8] .∈ (terms,))
        A = [A Bx_hat_Bx_dot Bx_hat_By_dot Bx_hat_Bz_dot]
        A = [A By_hat_Bx_dot By_hat_By_dot By_hat_Bz_dot]
        A = [A Bz_hat_Bx_dot Bz_hat_By_dot]
    end

    # add (3) eddy current terms - Bx_hat_Bx_dot, By_hat_By_dot, Bz_hat_Bz_dot
    if any([:eddy3,:e3] .∈ (terms,))
        A = [A Bx_hat_Bx_dot By_hat_By_dot Bz_hat_Bz_dot]
    end

    # add (3) derivative terms - Bx_dot, By_dot, Bz_dot
    if any([:fdm,:f,:d,:fdm3,:f3,:d3] .∈ (terms,))
        A = [A Bx_dot By_dot Bz_dot]
    end

    # add (1) bias term
    if any([:bias,:b] .∈ (terms,))
        A = [A ones(eltype(Bt),size(Bt))]
    end

    if return_B
        B_dot = [Bx_dot By_dot Bz_dot]
        return (A, Bt, B_dot)
    else
        return (A)
    end
end # function create_TL_A
*/

}

bool TollesLawson::fdm(std::vector<double> &dif, const std::vector<double> &x, const TollesLawson::FDMscheme scheme)
{
    //N = length(x)
		int N=x.size();
//		std::vector<double> dif(N);
		dif.resize(N);

		if(scheme==BACKWARD && N>1)
		{
			dif[0]=x[1]-x[0];
			for(int i=1;i<N;i++)
			{
				dif[i]=x[i]-x[i-1];
			}
			return true;
		}
		else if(scheme==FORWARD && N>1)
		{
			dif[N-1]=x[N-1]-x[N-2];
			for(int i=0;i<N-1;i++)
			{
				dif[i]=x[i+1]-x[i];
			}
			return true;
		}
		else if(scheme==CENTRAL&&N>2)
		{
			dif[0]=x[1]-x[0];
			dif[N-1]=x[N]-x[N-1];
			for(int i=1;i<N-1;i++)
			{
				dif[i]=0.5*(x[i+1]-x[i-1]);
			}
			return true;
		}
		else if(scheme==BACKWARD2&&N>3)
		{
			dif[0]=x[1]-x[0];
			dif[1]=x[2]-x[1];
			for(int i=2;i<N;i++)
			{
				dif[i]=0.5*(3.0*x[i]-4.0*x[i-1]+x[i-2]);
			}
			return true;
		}
		else if(scheme==FORWARD2 && N>3)
		{
			dif[N-1]=x[N-1]-x[N-2];
			dif[N-2]=x[N-2]-x[N-3];
			for(int i=0;i<N-2;i++)
			{
				dif[i]=0.5*(-x[i+2]+4.0*x[i+1]-3.0*x[i]);
			}
			return true;
		}
		else if(scheme==FOURTH && N>4)
		{
			dif[0]=0;
			dif[1]=0;
			dif[N-1]=0;
			dif[N-2]=0;
			for(int i=2;i<N-2;i++)
			{
				dif[i]=(x[i-2]-4.0*x[i-1]+6.0*x[i]-4.0*x[i+1]+x[i+2])/16.0;
			}
			return true;
		}
		else
		{
			std::fill(dif.begin(), dif.end(), 0.0);
			return false;
		}

    //if (scheme == :backward) & (N > 1)
    //    dif_1   =  x[2]       - x[1]
    //    dif_end =  x[end]     - x[end-1]
    //    dif_mid = (x[2:end-1] - x[1:end-2])
    //elseif (scheme == :forward) & (N > 1)
    //    dif_1   =  x[2]       - x[1]
    //    dif_end =  x[end]     - x[end-1]
    //    dif_mid = (x[3:end] - x[2:end-1])
    //elseif (scheme in [:central,:central2]) & (N > 2)
    //    dif_1   =  x[2]     - x[1]
    //    dif_end =  x[end]   - x[end-1]
    //    dif_mid = (x[3:end] - x[1:end-2]) ./ 2
    //elseif (scheme == :backward2) & (N > 3)
    //    dif_1   = x[2:3]     - x[1:2]
    //    dif_end = (3*x[end]     - 4*x[end-1]   + x[end-2]    ) ./ 2
    //    dif_mid = (3*x[3:end-1] - 4*x[2:end-2] + x[1:end-3]  ) ./ 2
    //elseif (scheme == :forward2) & (N > 3)
    //    dif_1   = (-x[3]        + 4*x[2]       - 3*x[1]      ) ./ 2
    //    dif_end = x[end-1:end] - x[end-2:end-1]
    //    dif_mid = (-x[4:end]    + 4*x[3:end-1] - 3*x[2:end-2]) ./ 2
    //elseif (scheme in [:fourth,:central4]) & (N > 4)
    //    dif_1   = zeros(eltype(x),2)
    //    dif_end = zeros(eltype(x),2)
    //    dif_mid = (   x[1:end-4] + 
    //               -4*x[2:end-3] + 
    //                6*x[3:end-2] + 
    //               -4*x[4:end-1] + 
    //                  x[5:end  ] ) ./ 16 # divided by dx^4
    //else
    //    return zero(x)
    //end

    //dif = [dif_1; dif_mid; dif_end]

    //return (dif)
}

} // namespace magnav

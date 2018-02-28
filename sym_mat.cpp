/* Author: Kushagra Chandak */

/* 
   SymMat class for storing symmetric matrices as 1-D array and manipulating data thereafter. 
   Data storage and access defined.
   Also, matrix multiplication and addtion defined with operators *, + and =
 */

/* Include files */
#include <iostream>
#include <assert.h>
#include <Eigen/Dense>
#include <ctype.h>

//template <class T, unsigned int N=3>
//Eigen::MatrixXd operator* (const SymMat<T,N>& lhs, const SymMat<T,N>& rhs);

template <class T, unsigned int N=3>
class SymMat
{
        public:
                SymMat(Eigen::MatrixXd M);

                /* Creates a 1-D array to store symmetric matrices */
                void setIndex(Eigen::MatrixXd M);

                /* Get index of the 1-D array */
                unsigned int getIndex(unsigned int i, unsigned int j);

                /* Only "()" operator is used to access the 1-D array using the 2-D matrix. Alternately, "[]" could also have been used */
                inline T& operator()(unsigned int i, unsigned int j)
                { return Array[getIndex(i, j)]; }

                inline SymMat<T, N>& operator=(const SymMat& rhs) 
                {
                        for(unsigned int i=0; i<size; ++i) Array[i] = rhs.Array[i];
                        return *this;
                }

                /* Addition of 2 symmetric matrices */
                friend SymMat<T, N> operator+(const SymMat<T,N>& lhs, const SymMat& rhs) 
                {
                        Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(N,N);
                        SymMat<double,3> add(temp);
                        for(unsigned int i=0; i<N*(N+1)/2; ++i) add.Array[i] = lhs.Array[i] + rhs.Array[i];
                        return add;
                }

                /* Addtion of a symmetric matrix and a generic matrix */
                friend Eigen::MatrixXd operator+(SymMat<T,N>& lhs, Eigen::MatrixXd& rhs)
                {
                        Eigen::MatrixXd add = Eigen::MatrixXd::Zero(N,N);
                        for(unsigned int i=0;i<N;i++)
                                for(unsigned int j=0;j<N;j++)
                                        add(i,j) = rhs(i,j)+lhs(i,j);
                        return add;
                }

                /* Multiplication of 2 matrices (Sym*Sym or Sym*Generic) */  
                template <class U>
                        friend Eigen::MatrixXd operator*(SymMat<T,N>& lhs, U& rhs)
                        {
                                Eigen::MatrixXd mult = Eigen::MatrixXd::Zero(N,N);
                                for(unsigned int i=0;i<N;i++)
                                {
                                        for(unsigned int j=0;j<N;j++)
                                        {
                                                for(unsigned int k=0;k<N;k++) mult(i,j) += lhs(i,k)*rhs(k,j);
                                        }
                                }
                                return mult;
                        }

        private:
                double* Array;

                /* size of the symmetric matrix (only the upper triangular and diagonal is stored)*/
                unsigned int size;
};

//template <class T, unsigned int N>
//Eigen::MatrixXd operator *(const SymMat<T,N>& lhs, const SymMat<T,N>& rhs)

template <class T, unsigned int N>
SymMat<T,N>::SymMat(Eigen::MatrixXd M)
{
        size = N*(N+1)/2;
        Array = new T[size];
        setIndex(M);
}

template <class T, unsigned int N>
void SymMat<T,N>::setIndex(Eigen::MatrixXd M)
{
        unsigned int index;
        for(unsigned int i=0;i<N;i++)
        {
                for(unsigned int j=0;j<N;j++)
                {
                        index = getIndex(i,j);
                        Array[index] = M(i,j);
                }
        }
}                       

template <class T, unsigned int N>
unsigned int SymMat<T,N>::getIndex(unsigned int row, unsigned int col)
{
        unsigned int index;
        if(row>col)
        {
                int temp = row;
                row = col;
                col = temp;
        }
        index = row*N - (row+1)*row/2 + col;
        return index;
}

bool isSymmetric(Eigen::MatrixXd M, unsigned int N=3)
{
        Eigen::MatrixXd MT = M.transpose();
        for(unsigned int i=0;i<N;i++)
        {
                for(unsigned int j=0;j<N;j++)
                        if(M(i,j) != MT(i,j)) return false;
        }
        return true;
}

void print(Eigen::MatrixXd M,unsigned int N=3)
{
        for(unsigned int i=0;i<N;++i)
        {
                for(unsigned int j=0;j<N;++j) std::cout << M(i,j) << " ";
                std::cout << std::endl;
        }
}

int main()
{
        unsigned int N=3;
        assert (N==3);
        double n;
        /* Taking 1st matrix as input */
        unsigned int count = 0;
        Eigen::MatrixXd M1(N,N);
        int i=0,j=0;
        std::cout << "Enter 1st (3x3) matrix (only numbers): " << std::endl;
        for(unsigned int i=0;i<N;i++)
        {
                for(unsigned int j=0;j<N;j++)
                        std::cin >> M1(i,j);
        }
        std::cout << std::endl;

        /* Taking 2nd matrix as input */
        Eigen::MatrixXd M2(N,N);
        std::cout << "Enter 2nd (3x3) matrix (only numbers): " << std::endl;
        for(unsigned int i=0;i<N;i++)
        {
                for(unsigned int j=0;j<N;j++)
                        std::cin >> M2(i,j);
        }
        std::cout << std::endl;
       

        /* Checking whether matrices are symmetric */
        bool symm1, symm2;
        symm1 = isSymmetric(M1,N);
        symm2 = isSymmetric(M2,N);
        SymMat<double,3> S1(M1);
        SymMat<double,3> S2(M2);

        /* Initializing result matrices */
        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
        Eigen::MatrixXd M_add(N,N);
        Eigen::MatrixXd M_mult(N,N);
        Eigen::MatrixXd M_symm_mult(N,N);
        SymMat<double,3> S(M);

       /*  Adding and multiplying the matrices */
        if(symm1 && symm2)
        {
                S=S1+S2;
                M_symm_mult=S1*S2;
                std::cout << "Addition of the 2 matrices is: " << std::endl;
                for(int i=0;i<N;i++)
                {
                        for(int j=0;j<N;j++)
                                std::cout << S(i,j) << " ";
                        std::cout << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Multiplication of 2 matrices is: " << std::endl;
                print(M_symm_mult,N);
        }
        else if(symm1 && !symm2)
        {
                M_add=S1+M2;
                M_mult=S1*M2;
                std::cout << "Addition of the 2 matrices is: " << std::endl;
                print(M_add,N);
                std::cout << std::endl;
                std::cout << "Multiplication of the 2 matrices is: " << std::endl;
                print(M_mult,N);
        }
        else if(!symm1 && symm2)
        {
                M_add=S2+M1;
                M_mult=S2*M1;
                std::cout << "Addition of the 2 matrices is: " << std::endl;
                print(M_add,N);
                std::cout << std::endl;
                std::cout << "Multiplication of the 2 matrices is: " << std::endl;
                print(M_mult,N);
        }
        else if(!symm1 && !symm2)
        {
                M_add=M1+M2;
                M_mult = M1*M2;
                std::cout << "Addition of the 2 matrices is: " << std::endl;
                print(M_add,N);
                std::cout << std::endl;
                std::cout << "Multiplication of the 2 matrices is: " << std::endl;
                print(M_mult,N);
        }

        return 0;
}


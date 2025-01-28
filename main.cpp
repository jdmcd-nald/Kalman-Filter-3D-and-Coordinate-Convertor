//
//  main.cpp
//  CSV_ECEF
//
//  Created by Jared on 12/26/24.
//

#include <iostream>
#include <fstream>
#include <filesystem>
#include <mdspan>

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <utility>
#include <string>

#include <array>
#include <iostream>
#include <any>

#include <cassert>

template <typename T, int N, int M>
int length(T(&)[N][M])
{
    return N;
}

class Matrix2D
{
private:
    // You might be tempted to make m_arr a reference to an ArrayFlat2d,
    // but this makes the view non-copy-assignable since references can't be reseated.
    // Using std::reference_wrapper gives us reference semantics and copy assignability.
    std::vector<long double> my2DArray; // 3x3 2D array
    int m_row;
    int m_column;

public:
    Matrix2D(int row,int col,std::vector<long double> vec) : m_row(row),m_column(col),my2DArray(vec) {}

    // Get element via single subscript (using operator[])
    //T& operator[](int i) { return m_arr[i]; }
    //const T& operator[](int i) const { return m_arr[i]; }

    // Get element via 2d subscript (using operator(), since operator[] doesn't support multiple dimensions prior to C++23)
    //T& operator[](int row, int col) { return m_arr[static_cast<std::size_t>(row * cols() + col)]; }
    //const T& operator[](int row, int col) const { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }

    
    
    // in C++23, you can uncomment these since multidimensional operator[] is supported
//    T& operator[](int row, int col) { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }
//    const T& operator[](int row, int col) const { return m_arr.get()[static_cast<std::size_t>(row * cols() + col)]; }
    
    void replace_data(Matrix2D& arr)
    {
        assert((*this).rows() == arr.rows());
        assert(       (*this).cols() == arr.cols());
        
        for(int k = 0;k<(*this).rows()*(*this).cols();k++)
        {
            ((*this).data())[k] = ((arr).data())[k];
        }
    }
    
    void scale(long double scaler)
    {
        
        std::vector<long double> output_vector;
        
        for(int k = 0;k<(*this).rows();k++)
        {
            for(int j = 0;j<(*this).cols();j++)
            {
                long double test_2 = scaler *((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_2;
            }
        }
    }
    
    void addition(Matrix2D& arr)
    {
        assert((*this).cols() == arr.cols());
        assert(       (*this).rows() == arr.rows());
        
        std::vector<long double> output_vector;
        
        for(int k = 0;k<(*this).rows();k++)
        {
            for(int j = 0;j<(*this).cols();j++)
            {
                long double test_1 = (arr.data())[k*(arr.cols())+j];
                long double test_2 = ((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_1+test_2;
            }
        }
    }
    
    void subtraction(Matrix2D& arr)
    {
        assert((*this).rows() == arr.rows());
        assert(       (*this).cols() == arr.cols());
        
        std::vector<long double> output_vector;
        
        for(int k = 0;k<(*this).rows();k++)
        {
            for(int j = 0;j<(*this).cols();j++)
            {
                long double test_1 = (arr.data())[k*(arr.cols())+j];
                long double test_2 = ((*this).data())[k*((*this).cols())+j] ;
                //need to replace with output_vector
                ((*this).data())[k*((*this).cols())+j] = test_2-test_1;
            }
        }
    }
    
    Matrix2D multiply(Matrix2D& arr)
    {
        
        assert((*this).cols() == arr.rows());
        
        std::vector<long double> output_vector;
        output_vector.reserve((*this).rows()*arr.cols());
        long double holder = 0;
        for(int j = 0;j<(*this).rows();j++)
        {
            for(int k = 0;k<(arr).cols();k++)
            {
                holder = 0;
                for(int i = 0;i<(*this).cols();i++)
                {
                    long double left = ((*this).data())[j*((*this).cols())+i];
                    long double right = (arr.data())[i*(arr.cols())+k];
                    holder += left*right;
                }
                output_vector.push_back(holder);
            }
        }
        Matrix2D output((*this).rows(),arr.cols(),output_vector);
        return output;
    }
    
    Matrix2D output_identity(int size)
    {
        std::vector<long double> output_vector;
        output_vector.reserve(size*size);
        for(int j = 0;j<size;j++)
        {
            for(int k = 0;k<size;k++)
            {
                if(j==k)
                {
                    output_vector.push_back(1);
                }
                else
                {
                    output_vector.push_back(0);
                }
            }
        }
        Matrix2D output(size,size,output_vector);
        return output;
    }
    
    Matrix2D transpose()
    {
        std::vector<long double> output_vector;
        output_vector.reserve((*this).rows()*(*this).cols());
        for(int j = 0;j<(*this).cols();j++)
        {
            for(int k = 0;k<(*this).rows();k++)
            {
                output_vector.push_back(((*this).data())[k*((*this).cols())+j]);
            }
        }
        Matrix2D output((*this).cols(),(*this).rows(),output_vector);
        return output;
    }
    
    
    
    Matrix2D cofactor()
    {
        assert((*this).rows() == (*this).cols());
        
        //std::vector<Matrix2D> cofactor_segment;
        //cofactor_segment.reserve((*this).rows() * (*this).cols());
        std::vector<long double> cofactor_matrix_vector;
        
        if(((*this).rows()) == 2)
        {
            cofactor_matrix_vector.push_back(((*this).data())[3]);
            cofactor_matrix_vector.push_back(-((*this).data())[1]);
            cofactor_matrix_vector.push_back(-((*this).data())[2]);
            cofactor_matrix_vector.push_back(((*this).data())[0]);
            Matrix2D cofactor_matrix((*this).rows(),(*this).cols(),cofactor_matrix_vector);
            return cofactor_matrix;
        }
        
        for(int j = 0;j<(*this).cols();j++)
        {
            for(int k = 0;k<(*this).rows();k++)
            {
                std::vector<long double> cofactor_segment_vector;
                cofactor_segment_vector.reserve((((*this).rows())-1) * (((*this).cols())-1));
                
                for(int i = 0;i<(*this).rows();i++)
                {
                    for(int o = 0;o<(*this).cols();o++)
                    {
                        if(i!=k && o!=j)
                        {
                            cofactor_segment_vector.push_back(((*this).data())[i*((*this).cols())+o]);
                        }
                    }
                }
                
                int u = -1;
                if(k%2 == j%2)
                {
                    u=1;
                }
                
                Matrix2D temp(((*this).rows())-1,((*this).cols())-1,cofactor_segment_vector);
                
                cofactor_matrix_vector.push_back(u*det(temp));
                
                //cofactor_segment.push_back(temp);
            }
        }
        
        Matrix2D cofactor_matrix((*this).rows(),(*this).cols(),cofactor_matrix_vector);
        return cofactor_matrix;
    }
    
    double long det(Matrix2D arr)
    {
        std::vector<long double> determinate;
        if((arr).rows()==2){
            return(det_2D(arr));
        }
        
        for(int j = 0;j<arr.cols();j++)
        {
            int k =0;
            {
                std::vector<long double> det_segment_vector;
                det_segment_vector.reserve((((arr).rows())-1) * (((arr).cols())-1));
                
                for(int i = 0;i<arr.rows();i++)
                {
                    for(int o = 0;o<arr.cols();o++)
                    {
                        if(i!=k && o!=j)
                        {
                            det_segment_vector.push_back(((arr).data())[i*((arr).cols())+o]);
                        }
                    }
                }
                
                int u = -1;
                if(k%2 == j%2){ u=1; }
                
                Matrix2D temp((((arr).rows())-1) ,(((arr).cols())-1),det_segment_vector);
                
                determinate.push_back(u*(((arr).data())[k*((arr).cols())+j])*det(temp));
                
                //cofactor_segment.push_back(temp);
            }
        }
        long double sum = 0;
        for(int i =0;i<determinate.size();i++)
        {
            sum += determinate[i];
        }
        return sum;
    }
    
    Matrix2D inverse(Matrix2D& arr)
    {
        Matrix2D output = ((arr.cofactor()).transpose());
        long double test = det(arr);
        output.scale(1.0/det(arr));
        return output;
    }
    
    long double det_2D(Matrix2D& arr)
    {
        assert(arr.rows() == arr.cols());
        return (arr.data())[0]*(arr.data())[3]-(arr.data())[1]*(arr.data())[2];
    }
    
    int rows() const { return m_row; }
    int cols() const { return m_column; }
    std::vector<long double>& data() { return my2DArray; }
};


class KalmanFilter3D
{
    long double dt = .10;
    long double process_noise_std;
    long double measurement_noise_std;
    //Matrix2D *A;
    
    Matrix2D H{3,6,{
        1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0}};
    
    Matrix2D P { H.output_identity(6)};
    Matrix2D P_result { H.output_identity(6)};
    
    Matrix2D A{ H.output_identity(6)};
    long double qq;
    Matrix2D Q{ H.output_identity(6)};
    
    Matrix2D R{ H.output_identity(3)};
    
    Matrix2D x{6,1,{0,0,0,0,0,0}};
    
public:
    KalmanFilter3D(long double deltat, long double p_noise_std, long double m_noise_std)
    : dt{deltat},process_noise_std{p_noise_std},measurement_noise_std{m_noise_std}
    {
        //Matrix2D *A;
        Matrix2D A_Matrix {6,6,
            {1, 0, 0, dt,0,0,
                0, 1, 0, 0,dt,0,
                0, 0, 1, 0,0,dt,
                0, 0, 0, 1,0,0,
                0, 0, 0, 0,1,0,
                0, 0, 0, 0,0,1}};
        
        A.replace_data(A_Matrix);
        
        qq = process_noise_std * process_noise_std;
        
        Matrix2D Q_Matrix (6,6,
            {  qq, 0, 0, 0, 0, 0,
                0, qq, 0, 0, 0, 0,
                0, 0, qq, 0, 0, 0,
                0, 0, 0, qq, 0, 0,
                0, 0, 0, 0, qq, 0,
                0, 0, 0, 0, 0, qq
            });
        
        Q.replace_data(Q_Matrix);
        
        long double r = measurement_noise_std * measurement_noise_std;
        
        Matrix2D R_Matrix(3,3, {
                r, 0, 0,
                0, r, 0,
                0, 0, r});
        
        R.replace_data(R_Matrix);
    }
    
    void predict(long double d_t)
    {
        
        Matrix2D A_Matrix {6,6,
            {1, 0, 0, d_t,0,0,
                0, 1, 0, 0,d_t,0,
                0, 0, 1, 0,0,d_t,
                0, 0, 0, 1,0,0,
                0, 0, 0, 0,1,0,
                0, 0, 0, 0,0,1}};
        
        A.replace_data(A_Matrix);
        
        Matrix2D x_result = A.multiply(x);
        x.replace_data(x_result);
        Matrix2D AP = A.multiply(P);
        Matrix2D At = A.transpose();
        Matrix2D APAt = AP.multiply(At);
        APAt.addition(Q);
        P.replace_data(APAt);
    }
    
    void predict()
    {
        Matrix2D x_result = A.multiply(x);
        x.replace_data(x_result);
        Matrix2D AP = A.multiply(P);
        Matrix2D At = A.transpose();
        Matrix2D APAt = AP.multiply(At);
        APAt.addition(Q);
        P.replace_data(APAt);
    }
    
    void update(Matrix2D z)
    {
        auto HP = H.multiply(P);
        
        auto Ht = H.transpose();
        auto HPHt = HP.multiply(Ht);
        
        HPHt.addition(R);
        auto S = HPHt;
        
        auto S_inv = S.inverse(S);
        
        auto PHt = P.multiply(Ht);
        
        auto K = PHt.multiply(S_inv);
        
        /////
        
        auto y_Hx = H.multiply(x);
        auto y_sub = H.multiply(x);
        auto y = z;
        y.subtraction(y_sub);
        
        
        auto Ky = K.multiply(y);
        x.addition(Ky);
        
        auto KH = K.multiply(H);
        
        
        auto I_KH= P.output_identity(6);
        I_KH.subtraction(KH);
        auto newP = I_KH.multiply(P);
        P.replace_data(newP);
        
    }
    
    void current_state()
    {
        std::cout << "Estimated position: (" << x.data()[0] << ", " << x.data()[1] << ", " << x.data()[2] << ")\n";
    }
    
};


class Constants
{
private:
    const long double pi = 3.14159;
    const long double a = 6378137;
    const long double f = 1/298.257223563;
    const long double b = a*(1-f);
    const long double e = std::sqrtl((std::pow(a,2.)-std::pow(b,2.))/std::pow(a,2.));
    const long double ep = std::sqrtl((std::pow(a,2.)-std::pow(b,2.))/std::pow(b,2.));
public:
    // only get functions no sets and auto to make things shorter
    auto get_pi(){return pi;}
    auto get_a(){return a;}
    auto get_f(){return f;}
    auto get_b(){return b;}
    auto get_e(){return e;}
    auto get_ep(){return ep;}
};


struct Units
{
private:
    //angle conversion
    const long double ang_rad = 1;
    const long double ang_deg = 3.14159/180;
    //length conversion
    const long double km = 1000;
    const long double m = 1;
    const long double cm = .001;
    const long double mm = .0001;
public:
    
    double convert(std::string current_unit,std::string new_unit)
    {
        return (exclusion_table(current_unit).first)/(exclusion_table(new_unit).first);
    }
    
    std::pair<long double, std::string> exclusion_table(std::string unit)
    {
        if(unit=="kilometer") return std::make_pair(km,"length");
        if(unit=="meter") return std::make_pair(m,"length");
        if(unit=="centimeter") return std::make_pair(cm,"length");
        if(unit=="millimeter") return std::make_pair(mm,"length");
        if(unit=="radian") return std::make_pair(ang_rad,"angle");
        if(unit=="degree") return std::make_pair(ang_deg,"angle");
        std::cout<<"unknown unit conversion";
        return std::make_pair(1,"unk");
    }
};

//class to hold csv object
class CSV_Data
{
private:
    std::string m_file_name;
    size_t m_rows;
    size_t m_total_number_of_rows;
    //holds the data read from the csv
    std::vector<long double> data;
    
    //performs the actual splitting of the getline string for csv read
    std::vector<std::string> split(const std::string &input_string, std::string_view delim)
    {
        size_t start_pos = 0;
        size_t end_pos = 0;
        std::string sub_str;
        //vector which holds the split contents
        std::vector<std::string> s_vector;
        
        while ((end_pos = input_string.find(delim, start_pos)) < input_string.npos) {
            sub_str = input_string.substr(start_pos, end_pos - start_pos);
            start_pos = end_pos + delim.length();
            s_vector.push_back(sub_str);
        }
        //gets the final value after the last comma
        s_vector.push_back(input_string.substr(start_pos));
        return s_vector;
    }
    
public:
    CSV_Data(std::string name)
    : m_file_name{name}
    {
        std::string line;
        std::ifstream file;
        file.open(m_file_name);
        std::string delimiter = ",";
        size_t j=0;
        while(getline(file, line))
        {
            j++;
            std::vector<std::string> v = split(line, delimiter);
            m_rows = v.size();
            for (auto i : v)
            {
                //long double conversion used for percision
                //stold also removes \n and trailing spaces so need for
                //other checks to remove extra characters
                data.push_back(std::stold(i));
               // std::cout <<  std::stold(i) << ",";
            }
           // std::cout <<  std::endl;
        }
        m_total_number_of_rows = j; //easiest to calculate here and simplifies other functions by allowing reserve to be called
    }
    
    //get functions for private members
    std::vector<long double> get_data()
    {
        return data;
    }
    size_t get_rows()
    {
        return m_rows;
    }
    size_t get_max_rows()
    {
        return m_total_number_of_rows;
    }
    
};


class Base_Data //base class for LLA and ECEF data, reduce code bloat and duplication
{
public: // values were private initial but are now derived classes from Base_Data, so members
    // are protected now to stil be accessible by derived classes
    std::vector<long double> m_time;
    size_t m_number_of_rows;
    size_t m_initial_data_rows;
    
    std::vector <long double> x;
    std::vector <long double> y;
    std::vector <long double> z;
    std::vector <long double> x_vel;
    std::vector <long double> y_vel;
    std::vector <long double> z_vel;
    
    Units unit;
    
    //Constants can be looked to by anyone and it coontains no sets, so multiple owners is encoouraged to not waste memory
    std::shared_ptr<Constants> m_constants = std::make_shared<Constants>();

    //linear interpolation of data points
    long double lin_interp(const long double ukn_t,const std::vector<long double> &t_data,
                           const std::vector<long double> &y_data, const size_t cur_row, const size_t y_row)
    {
        long double y0 = y_data[3*(cur_row-1)+y_row];
        long double x = ukn_t;
        long double x0 = t_data[cur_row-1];
        long double y1 = y_data[cur_row*3+y_row];
        long double x1 = t_data[cur_row];
        return y0+(x-x0)*((y1-y0)/(x1-x0));
    }
    
    // function that reads data from a csv to pos and vel data vectors of the derived classes
    // data is imported as const as it wont change
    // pos_data and vel_data are imported as reference so that the values in the derived classes will be properly updated
    virtual void read_data_in(const std::vector<long double> &data , size_t number_of_rows,std::vector<long double> &pos_data, std::vector<long double> &vel_data)
    {
        m_number_of_rows=number_of_rows;
        m_initial_data_rows = data.size()/number_of_rows;
        
        
        
        KalmanFilter3D kf (.1,.1, .5);
        
        //if the data only contains time and pos data this is run
        if(m_initial_data_rows == 4)
        {
            m_time.reserve(number_of_rows);
            pos_data.reserve(m_number_of_rows);
            std::vector <std::vector<long double>*> pos_xyz_data = {&x,&y,&z};
            vel_data.reserve(m_number_of_rows);
            std::vector <std::vector<long double>*> pos_xyz_vel_data = {&x_vel,&y_vel,&z_vel};
            for(size_t i = 0;i<data.size();i++)
            {
                if(i%m_initial_data_rows==0)
                {
                    m_time.push_back(data[i]);
                }
                else if(i%m_initial_data_rows>=1 and i%m_initial_data_rows<=3)
                {
                    pos_data.push_back(data[i]);
                    
                    pos_xyz_data[(i%m_initial_data_rows)-1]->push_back(data[i]);
                    long double qq=(*pos_xyz_data[0])[0];
                    if(x[0] == pos_xyz_data[0]->at(0))
                    {
                        x[0]++;
                        x[0]--;
                    }
                    
                    //calculation for whiiich row is currently being indexed
                    //there are faster methods such as increment a variable once
                    //another exceeds m_initial_data_rows
                    int row = (int(i)-int(i%m_initial_data_rows))/m_initial_data_rows;
                    if(i>m_initial_data_rows)
                    {
                        //velocity calc
                        long double time_dif = (m_time[(row-1)]-m_time[row]);
                        long double pos_dif = (pos_data[(row-1)*3+((i%4)-1)]-pos_data[row*3+((i%4)-1)]);
                        vel_data.push_back(pos_dif/time_dif);
                        pos_xyz_vel_data[(i%m_initial_data_rows)-1]->push_back(pos_dif/time_dif);
                    }
                    else
                    {
                        //for the first position a veloocity of 0 is provided
                        //though this could later be interpolated from future data points
                        vel_data.push_back(0);
                        pos_xyz_vel_data[(i%m_initial_data_rows)-1]->push_back(0);
                    }
                }
                if(i%m_initial_data_rows == 3)
                {
                    if(i<4)
                    {
                        Matrix2D KalmanTest(3,1,
                            {x[0],
                            y[0],
                            z[0]});
                        kf.predict();
                        kf.update(KalmanTest);
                    }
                    else
                    {
                        Matrix2D KalmanTest(3,1,
                            {x[i/4],
                            y[i/4],
                            z[i/4]});
                        kf.predict();
                        kf.update(KalmanTest);
                        kf.current_state();
                    }
                    
                }
            }
        }
        
        //if data contains pos and vel data this is run
        if(m_initial_data_rows == 7)
        {
            m_time.reserve(number_of_rows);
            pos_data.reserve(number_of_rows);
            vel_data.reserve(number_of_rows);
            for(size_t i=0; i<data.size();i++)
            {
                if(i%7==0)
                {
                    m_time.push_back(data[i]);
                }
                else if(i%7>=1 and i%7<=3)
                {
                    pos_data.push_back(data[i]);
                }
                else if(i%7>=4 and i%7<=6)
                {
                    vel_data.push_back(data[i]);
                }
            }
        }
    }

    //function to interpolate properties at various times
    virtual std::vector<long double> interpolate_at_time(long double requested_time,const std::vector<long double> &input_data)
    {
        std::vector<long double> output;
        output.reserve(3);
        for(size_t i = 0;i<m_time.size();i++)
        {
            if(m_time[0]>requested_time) //allows interpolation before the start time though accuracy will be expectedly poor
            {
                std::cout << "veracity of values before inputted data should not be relied upon" << std::endl;
                for(size_t j=0;j<3;j++)
                {
                    output.push_back(lin_interp(requested_time,m_time,input_data,i+1,j));
                    std::cout << lin_interp(requested_time,m_time,input_data,i+1,j) << ",";
                }
                std::cout <<std::endl;
                return output;
            }
            if(m_time[i]>=requested_time)
            {
                for(size_t j=0;j<3;j++)
                {
                    output.push_back(lin_interp(requested_time,m_time,input_data,i,j));
                    std::cout << lin_interp(requested_time,m_time,input_data,i,j) << ",";
                }
                std::cout << std::fixed << std::setprecision(0) << requested_time << std::endl;
                return output;
            }
        }
        return output;
    }
public:
    //gets for rows of the initial input data
    //gets for total number of rows
    size_t get_rows(){return m_initial_data_rows;}
    size_t get_max_rows(){return m_number_of_rows;}
    
    //virtual functions for property interpolation
    //two are used for disambiguation as the velocity or position could
    //be calculated in different way with added data but one function is also appropriate here
    
    virtual std::vector<long double> pos_at_time(long double requested_time) = 0;
    virtual std::vector<long double> vel_at_time(long double requested_time) = 0;
};

class ECEF : virtual public Base_Data
{
    //xyz values for pos and vel
    //not in base class as technically lla vel gives purely angular changes
    //while xyz vel gives linear
    
    std::vector<long double>* x;
    std::vector<long double>* y;
    std::vector<long double>* z;
    std::string x_units = "meter";
    std::string y_units = "meter";
    std::string z_units = "meter";
    
    std::vector<long double> m_xyz_pos_data;
    std::vector<std::string> pos_units = {"meter","meter","meter"};
    std::vector<long double> m_xyz_vel_data;
    std::vector<std::string> vel_units = {"meters/s","meters/s","meters/s"};
public:
    ECEF(const std::vector<long double> &data , size_t number_of_rows)
    {
        //data read in
        Base_Data::read_data_in(data,number_of_rows, m_xyz_pos_data, m_xyz_vel_data);
    }
    
    std::vector<long double> pos_at_time(long double requested_time)
    {
        return Base_Data::interpolate_at_time(requested_time,m_xyz_pos_data);
    }
    std::vector<long double> vel_at_time(long double requested_time)
    {
        return Base_Data::interpolate_at_time(requested_time,m_xyz_vel_data);
    }
    
    std::vector<long double> Convert_to_LLA()
    {
        //conversion to LLA
        std::vector<long double> LLA_data;
        LLA_data.reserve(m_number_of_rows*4);
        
        //These gets and variable definiattiions are used
        //to simplify writing the calculatins and could be added to Base_Data
        //but they feel more appropriate close to  their usage
        //if more calculations were performed they should be moved to the base class
        //In a sipmler proogram these could alsoo be DEFINE constants
        const long double pi = m_constants->get_pi();
        const long double a = m_constants->get_a();
        const long double b = m_constants->get_b();
        const long double e = m_constants->get_e();
        const long double e2 = std::pow(e,2.);
        const long double ep = m_constants->get_ep();
        const long double ep2 = std::pow(ep,2.);
        
        for( size_t i =0;i< m_number_of_rows;i++)
        {
            //stores XYZ
            long double X = m_xyz_pos_data[i*3+0];
            long double Y = m_xyz_pos_data[i*3+1];
            long double Z = m_xyz_pos_data[i*3+2];
            
            long double lambda = (180./pi)*atan2(Y,X);
            long double p = std::sqrtl(pow(X,2)+pow(Y,2));
            
            //closed form solution
            long double theta = std::atan2((Z*a),(p*b));
            long double phi = (180./pi)*std::atan2((Z+ep2*b*std::pow(std::sinl(theta),3)),
                                   (p- e2*a*std::pow(std::cosl(theta),3)));// check
            long double N=a/std::sqrtl(1-e2*std::pow(std::sinl(phi),2.));
            long double h = (p/std::cosl((p/180)*phi))-N;
            //h can vary significantly by the percision of phi angle
            //and the resulting h valuue has variabillity +/- 1000 m from the real solution
            
            
            //iterative method
            long double h0=0.0000001;
            long double hi=LONG_MAX;
            long double phi0 = std::atan(Z/(p*(1-e2)));
            long double phii=LONG_MAX;
            int iii =0;
            //The percent difference between the i and i+1 values is used as the convergence criterea
            //A percent difference of.1% is used below
            while(100*(std::abs((h0-hi))/((h0+hi)/2))>.1 ||
                  100*(std::abs((phi0-phii))/((phi0+phii)/2))>.1)
            {
                if(iii>0){
                    h0=hi; phi0=phii;
                }
                N=a/std::sqrtl(1-e2*std::pow(std::sinl(phi0),2.));
                hi = (p/std::cosl(phi0)) - N;
                phii = atan(Z/(p*(1-e2*(N/(N+hi)))));
                iii++;
            }
            //iterative method shows much better convergence to actual h
            //when verified by the initial input data
            
            phi = (180./pi)*phii;
            h = hi;
            LLA_data.push_back(Base_Data::m_time[i]);
            LLA_data.push_back(phi);
            LLA_data.push_back(lambda);
            LLA_data.push_back(h);
            
            //print out to console
            std::cout << std::fixed
             << std::setprecision(std::numeric_limits<double>::max_digits10) << Base_Data::m_time[i]<<","<<phi<<","<<
            lambda<<"," << h << std::endl;
        }
        return LLA_data;
    }
    
};

class LLA : virtual public Base_Data
{
    //lattitude and longitude and altitude positions and
    //angular velocities and linear altitude velocities
    std::vector<long double> m_lla_pos_data;
    std::vector<long double> m_lla_vel_data;
    
    std::vector<long double>* lat;
    std::vector<long double>* lon;
    std::vector<long double>* alt;
    std::string lat_unit = "degree";
    std::string lon_unit = "degree";
    std::string alt_unit = "meter";
    
    std::vector<long double>* lat_rate;
    std::vector<long double>* lon_rate;
    std::vector<long double>* alt_rate;
    std::string lat_rate_unit = "degree";
    std::string lon_rate_unit = "degree";
    std::string alt_rate_unit = "meter";
    
    std::vector<long double> altalt;
    std::vector<std::string> pos_units = {"rad","rad","meter"};
    std::vector<std::string> vel_units = {"rads/s","rads/s","rads/s"};
    
    
public:
    LLA(const std::vector<long double> &data, size_t number_of_rows)
    {
        Base_Data::read_data_in(data, number_of_rows, m_lla_pos_data, m_lla_vel_data);
        lat = &(Base_Data::x);
        lon = &(Base_Data::y);
        alt = &(Base_Data::z);
        altalt = (Base_Data::x);
        long double bb = altalt[0];
        bb = alt->at(0);
        bb = (*alt)[0];
        long double r = unit.convert("meter", "millimeter");
        r = bb;
    }
    
    std::vector<long double> pos_at_time(long double requested_time)
    {
        return Base_Data::interpolate_at_time(requested_time,m_lla_pos_data);
    }
    std::vector<long double> vel_at_time(long double requested_time)
    {
        return Base_Data::interpolate_at_time(requested_time,m_lla_vel_data);
    }
    
    std::vector<long double> Convert_to_ECEF()
    {
        //conversion to ECEF
        std::vector<long double> ECEF_data;
        ECEF_data.reserve(m_number_of_rows*4);
        const long double pi = m_constants->get_pi();
        const long double a = m_constants->get_a();
        const long double a2 = std::pow(a,2.);
        const long double b = m_constants->get_b();
        const long double b2 = std::pow(b,2.);
        const long double e = m_constants->get_e();
        const long double e2 = std::pow(e,2.);
        long double phi ;
        long double lambda;
        long double h;
        long double N;
        for( size_t i =0;i< m_number_of_rows;i++)
        {
            //C++23 mdspan example. All data is currently stored a 1D vector
            //but mmdspan allows for multidimensional representations
            //has no impact on the actual math, but it can assist with visualizatioon
            auto ms2 = std::mdspan(m_lla_pos_data.data(), m_number_of_rows, 3);
            
            //std::cout << ms2[0,0] <<","<<ms2[0,1]<<","<<ms2[0,2]<<std::endl;
            //std::cout << ms2[1,0] <<","<<ms2[1,1]<<","<<ms2[1,2]<<std::endl;
            
            phi = (pi/180.)*(m_lla_pos_data[i*3+0]);
            phi = ((*lat)[i])*unit.convert(lon_unit,"radian");
            lambda = (pi/180.)*(m_lla_pos_data[i*3+1]);
            lambda = ((*lon)[i])*unit.convert(lat_unit,"radian");
            h = m_lla_pos_data[i*3+2];
            h = ((*alt)[i]);
            
            
            //equivalent to above
            phi = (pi/180.)*(ms2[i,0]);
            lambda = (pi/180.)*(ms2[i,1]);
            h = ms2[i,2];
            
            //calculation of phi
            N=a/std::sqrtl(1-e2*std::pow(std::sinl(phi),2.));
            ECEF_data.push_back(Base_Data::m_time[i]);
            ECEF_data.push_back((N+h)*cosl(phi)*cosl(lambda));
            ECEF_data.push_back((N+h)*cosl(phi)*sinl(lambda));
            ECEF_data.push_back(((b2/a2)*N+h)*std::sinl(phi));
            
            //velocities could also be calculated and appended to the vectors here but that would require recalculation
            //the data with added velocccities would then be read in as a 7 row data list.
            
            std::cout << Base_Data::m_time[i]<<","<<(N+h)*cosl(phi)*cosl(lambda)<<","<<
             (N+h)*cosl(phi)*sinl(lambda)<<"," <<
              ((b2/a2)*N+h)*std::sinl(phi)<< std::endl;
        }
        
        return ECEF_data;
    }
    
};

int main(int argc, const char * argv[]) {
    
    KalmanFilter3D kf (.1,.1, .5);
    
    {
        Matrix2D KalmanTest(3,1,
                            {-0.0161335, 0.284842, -0.29425});
        kf.predict();
        kf.update(KalmanTest);
        kf.current_state();
    }
    {
        Matrix2D KalmanTest(3,1,
                         {-0.249067, 0.220212, 0.159958});
        kf.predict();
        kf.update(KalmanTest);
        kf.current_state();
    }
    {
        Matrix2D KalmanTest(3,1,
                         {0.106394, -0.334886, -0.537161});
        kf.predict();
        kf.update(KalmanTest);
        kf.current_state();
    }
    {
        Matrix2D KalmanTest(3,1,
                         {0.438168, -0.501777, -0.762506});
        kf.predict();
        kf.update(KalmanTest);
        kf.current_state();
    }
    {
        Matrix2D KalmanTest(3,1,
                         {0.212342, 1.02412, 1.063343});
        kf.predict();
        kf.update(KalmanTest);
        kf.current_state();
    }
    
    
    
    Matrix2D cofactor_test_2(2,2,
        {1,2,
         3,4});
    
    cofactor_test_2.cofactor();
    
    Matrix2D cofactor_test_3(3,3,
        {1,2,3,
         4,5,6,
         7,8,9});
    
    cofactor_test_3.cofactor();
    
    Matrix2D cofactor_test_33(3,3,
        {1,2,3,
         2,5,4,
         3,3,2});
    
    cofactor_test_33.cofactor();
    
    
    Matrix2D cofactor_test_4(4,4,
        {1,2,3,4,
         2,4,1,1,
         3,4,1,2,
         4,2,1,3});
    
    cofactor_test_4.cofactor();
    
    
    cofactor_test_4.inverse(cofactor_test_4);
    
    
    Matrix2D test(3,3,
        {6,2,4,
        -1,4,3,
        -2,9,3});
    
    Matrix2D test1(3,1,
        {4,
        -2,
         1});
    
    Matrix2D testt(3,3,
        {2,3,2,
         2,3,1,
         2,1,0});
    
    Matrix2D testt1(3,2,
        {3,1,
         3,2,
         1,0});
    
    test.multiply(test1);
    testt.multiply(testt1);
    auto id = testt.output_identity(3);
    
    Matrix2D transpose_test(2,3,
        {1,2,3,
         4,5,6});
    
    transpose_test.transpose();
    
    Matrix2D transpose_test_2(3,3,
        {1,2,3,
        4,5,6,
        7,8,9});
    
    transpose_test_2.transpose();
    
    Matrix2D addition_test(2,2,
        {8,5,
        2,3});
    Matrix2D addition_test_2(2,2,
        {9,5,
        1,2});
    addition_test.addition(addition_test_2);
    
    Matrix2D cofactor_test(3,3,
        {11,12,13,
         21,22,23,
         31,32,33});
    
    std::string file_path = "LLA_Data.csv";
    if (std::filesystem::exists(file_path)) {
        std::cout << "File exists" << std::endl;
    } else {
        file_path = "/Users/user/Documents/CPP Projects/U/ECEF_Interview_Test/CSV_ECEF/CSV_ECEF/LLA_Data.csv";
        if (std::filesystem::exists(file_path)) {
                std::cout << "File exists." << std::endl;
        } else {
                std::cout << "File does not exist." << std::endl;
        }
    }
    
    CSV_Data read_data(file_path);

    LLA lla_data(read_data.get_data(),read_data.get_max_rows());
    
    ECEF ecef_data(lla_data.Convert_to_ECEF(),lla_data.get_max_rows());
    
    //velocities interpolated at requested Unix Times
    ecef_data.vel_at_time(1532334000);
    ecef_data.vel_at_time(1532335268);
    
    ecef_data.Convert_to_LLA(); // valiidation for ECEF transformations
    
    return 0;
}

#ifndef _TSP_position_h
#define _TSP_position_h


class TSP_position {

public: TSP_position();
	TSP_position(double x, double y);
	~TSP_position();

	void SetX(double x) {m_x=x;};
	void SetY(double y) {m_y=y;};
	double GetX() const {return m_x;};
	double GetY() const {return m_y;};

private: double m_x, m_y;		 

};

#endif

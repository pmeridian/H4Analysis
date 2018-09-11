#ifndef __FITPIX_RECO__
#define __FITPIX_RECO__

#include "interface/PluginBase.h"
#include "interface/FitpixTree.h"


using namespace std;

#define FITPIX_PIXELS_X 256
#define FITPIX_PIXELS_Y 256

class FitpixReco: public PluginBase
{
public:

  //---ctors---
  FitpixReco() {};
  
  //---dtor---
  ~FitpixReco() {};
  
  class FPHit {
  public:
    FPHit( int x, int y, float c=1 )
      {
	x_=x;
	y_=y;
	c_=c;
      }

    ~FPHit() {};

    inline int deltax( const FPHit& otherhit ) const
    {
      return x_-otherhit.x_;
    }
    inline int deltay( const FPHit& otherhit ) const
    {
      return y_-otherhit.y_;
    }
    inline bool isAdjacent( const FPHit& otherhit ) const
    {
      return ( (fabs(deltax(otherhit))<=1) && (fabs(deltay(otherhit))<=1) );
    }
    inline float deltaR( const FPHit& otherhit ) const
    {
      return sqrt( (float)deltax(otherhit)*deltax(otherhit) + (float)deltay(otherhit)*deltay(otherhit) );
    }

    int x_;
    int y_;
    float c_;
  };

  class FPCluster {
  public:
    
    FPCluster( const std::vector<FPHit>& hits )
      {
	hits_=hits;
      }
    FPCluster()
      {
	hits_.clear();
      }
    
    ~FPCluster() {};
    
    void add_hit( const FPHit& hit )
    {
      this->hits_.push_back(hit);
    }

    float x() const
    {
      float num(0.), denom(0.);

      for( unsigned i=0; i<hits_.size(); ++i ) {

	num   += hits_[i].c_ * hits_[i].x_;
	denom += hits_[i].c_;

      }

      return num/denom;

    }

    float y() const
    {
      float num(0.), denom(0.);

      for( unsigned i=0; i<hits_.size(); ++i ) {

	num   += hits_[i].c_ * hits_[i].y_;
	denom += hits_[i].c_;

      }

      return num/denom;

    }

    float charge() const
    {
      float charge(0.);

      for( unsigned i=0; i<hits_.size(); ++i ) {
	charge  += hits_[i].c_ ;
      }

      return charge;
    }

    bool isAdjacent( const FPHit& hit ) const
    {
      bool isAdj = false;

      for( unsigned i=0; i<this->hits_.size(); ++i ) {

	if( hit.isAdjacent( this->hits_[i] ) ) {
  
	  isAdj = true;
	  break;

	} // if

      } // for 

      return isAdj;
    }

    inline int nhits() const 
    { 
      return hits_.size(); 
    }

    int width() const
    {
      return std::max( this->widthx(), this->widthy() );
    }

    int widthx() const
    {
      int max_x = -1;
      int min_x = 9999;

      for( unsigned i=0; i<this->hits().size(); ++i ) {

	if( hits()[i].x_ > max_x ) max_x = hits()[i].x_;
	if( hits()[i].x_ < min_x ) min_x = hits()[i].x_;

      }

      return (max_x-min_x);

    }

    int widthy() const
    {
      int max_y = -1;
      int min_y = 9999;

      for( unsigned i=0; i<this->hits().size(); ++i ) {

	if( hits()[i].y_ > max_y ) max_y = hits()[i].y_;
	if( hits()[i].y_ < min_y ) min_y = hits()[i].y_;

      }

      return (max_y-min_y);
    }

    float rms() const
    {
      float xm = this->x();
      float ym = this->y();

      float var_x(0.);
      float var_y(0.);

      float sumW(0.);

      for( unsigned i=0; i<hits_.size(); ++i ) {

	sumW += hits_[i].c_;

	var_x += ( hits_[i].x_ - xm )*( hits_[i].x_ - xm )*hits_[i].c_; 
	var_y += ( hits_[i].y_ - ym )*( hits_[i].y_ - ym )*hits_[i].c_; 

      }

      var_x /= sumW;
      var_y /= sumW;

      float rms = std::sqrt( var_x + var_y );
  
      return rms;

    }

    float dist( const FPCluster& other_cluster ) const
    {
      float x1 = this->x();
      float y1 = this->y();

      float x2 = other_cluster.x();
      float y2 = other_cluster.y();

      return std::sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    }

    std::vector<FPHit> hits() const { return hits_; };
    void set_hits( const std::vector<FPHit>& hits ) { hits_ = hits;};

  private:

    std::vector<FPHit> hits_;

  };

    //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool ProcessEvent(const H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    
private:
    long int        boardId_;
    std::vector<FPHit> hits_;
    std::vector<FPCluster> clusters_;
    FitpixTree     fitpixTree_;
};

DEFINE_PLUGIN(FitpixReco);

#endif

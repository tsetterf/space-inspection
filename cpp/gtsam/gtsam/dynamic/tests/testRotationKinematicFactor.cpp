/**
 * @file   testRotationKinematicFactor.cpp
 * @brief  Unit tests for RotationKinematicFactor class
 * @author Tim Setterfield
 **/

#include <Eigen/Core>
#include <gtsam/dynamic/RotationKinematicFactor.h>
#include <gtsam/base/Testable.h>
#include <CppUnitLite/TestHarness.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Unit3.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/base/numericalDerivative.h>

using namespace std;
using namespace gtsam;

GTSAM_CONCEPT_TESTABLE_INST(RotationKinematicFactor)

static const Vector3 phiRBitoGi( (Vector(3) <<0.1952901,0.5282155,-0.0529983).finished() );
static const Vector3 phiRBjtoGj( (Vector(3) <<-0.0000000,0.0000000,0.6900000).finished() );
static const Vector3 phiRBitoW( (Vector(3) <<0.2000000,0.0000000,-0.0000000).finished() );
static const Vector3 phiRBjtoW( (Vector(3) <<-0.0000000,0.0000000,0.7900000).finished() );
static const Point3  tGitoBi_Gi( (Vector(3) <<0.6128071,0.0000000,-0.5055333).finished() );
static const Point3  tGjtoBj_Gj( (Vector(3) <<0.3121633,0.5600284,-0.0200000).finished() );
static const Point3  tGtoB_G( (Vector(3) <<-0.2500000,0.0000000,0.0000000).finished() );
static const Point3  tWtoBi_W( (Vector(3) <<1.0000000,0.0000000,0.0000000).finished() );
static const Point3  tWtoBj_W( (Vector(3) <<0.5038453,0.6103533,0.0000000).finished() );
static const Point3 delT( (Vector(3) <<0.0004000,-0.0030000,0.0200000).finished() );

//******************************************************************************
TEST(RotationKinematicFactor, Concept) {
  BOOST_CONCEPT_ASSERT((IsTestable<RotationKinematicFactor>));
}

//******************************************************************************
TEST(RotationKinematicFactor, equals) {
  Symbol p1('p',1), p2('p',2), p3('p',3), p4('p',4), t1('t',1);
  SharedNoiseModel noiseMod = noiseModel::Diagonal::Sigmas(
      (Vector(3) << 0.01, 0.01, 0.01).finished() );
  RotationKinematicFactor rkf1(p1, p2, p3, p4, t1, noiseMod);
  RotationKinematicFactor rkf2(p1, p2, p3, p4, t1, noiseMod);
  EXPECT(rkf1.equals(rkf2));
  RotationKinematicFactor rkf3;
  EXPECT(!rkf1.equals(rkf3));
}

//******************************************************************************
TEST(RotationKinematicFactor , Constructor) {
  Symbol p1('p',1), p2('p',2), p3('p',3), p4('p',4), t1('t',1);
  SharedNoiseModel noiseMod = noiseModel::Diagonal::Sigmas(
      (Vector(3) << 0.01, 0.01, 0.01).finished() );
  RotationKinematicFactor rkf1(p1, p2, p3, p4, t1, noiseMod);
  RotationKinematicFactor rkf2;
}

//******************************************************************************
TEST(RotationKinematicFactor , evaluateError) {
  Symbol p1('p',1), p2('p',2), p3('p',3), p4('p',4), t1('t',1);
  SharedNoiseModel noiseMod = noiseModel::Diagonal::Sigmas(
      (Vector(3) << 0.01, 0.01, 0.01).finished() );
  RotationKinematicFactor rkf1(p1, p2, p3, p4, t1, noiseMod);
  Rot3 RBitoW(Rot3::Rodrigues(phiRBitoW));
  Rot3 RBjtoW(Rot3::Rodrigues(phiRBjtoW));
  Rot3 RBitoGi(Rot3::Rodrigues(phiRBitoGi));
  Rot3 RBjtoGj(Rot3::Rodrigues(phiRBjtoGj));
  Pose3 PWtoBi(RBitoW,tWtoBi_W);
  Pose3 PWtoBj(RBjtoW,tWtoBj_W);
  Pose3 PGitoBi(RBitoGi,tGitoBi_Gi);
  Pose3 PGjtoBj(RBjtoGj,tGjtoBj_Gj);
  Matrix H1, H2, H3, H4, H5;

  Vector3 err = rkf1.evaluateError(PWtoBi, PGitoBi, PWtoBj, PGjtoBj,
		  tGtoB_G, H1, H2, H3, H4, H5);

  // Check that error matches the error artificially inserted
  EXPECT(assert_equal(delT,Point3(err),1e-6));

}

/* ************************************************************************* */
Vector3 rkfErr(RotationKinematicFactor& rkf, const Pose3& P1,
  const Pose3& P2, const Pose3& P3, const Pose3& P4, const Point3& t) {
  return rkf.evaluateError(P1, P2, P3, P4, t, boost::none, boost::none,
    boost::none, boost::none,boost::none);
}
TEST(RotationKinematicFactor , jacobians) {

  Symbol p1('p',1), p2('p',2), p3('p',3), p4('p',4), t1('t',1);
  SharedNoiseModel noiseMod = noiseModel::Diagonal::Sigmas(
      (Vector(3) << 0.01, 0.01, 0.01).finished() );
  RotationKinematicFactor rkf1(p1, p2, p3, p4, t1, noiseMod);
  Rot3 RBitoW(Rot3::Rodrigues(phiRBitoW));
  Rot3 RBjtoW(Rot3::Rodrigues(phiRBjtoW));
  Rot3 RBitoGi(Rot3::Rodrigues(phiRBitoGi));
  Rot3 RBjtoGj(Rot3::Rodrigues(phiRBjtoGj));
  Pose3 PWtoBi(RBitoW,tWtoBi_W);
  Pose3 PWtoBj(RBjtoW,tWtoBj_W);
  Pose3 PGitoBi(RBitoGi,tGitoBi_Gi);
  Pose3 PGjtoBj(RBjtoGj,tGjtoBj_Gj);
  Matrix H1, H2, H3, H4, H5;

  Vector3 err = rkf1.evaluateError(PWtoBi, PGitoBi, PWtoBj, PGjtoBj,
          tGtoB_G, H1, H2, H3, H4, H5);

  // Constant for numerical differentiation
  double delta = 1e-5;

  // Get expected jacobians using numerical differentiation
  Matrix H1e = numericalDerivative11<Vector3, Pose3>(
      boost::bind(rkfErr, rkf1, _1, PGitoBi, PWtoBj, PGjtoBj, tGtoB_G),
        PWtoBi, delta);
  Matrix H2e = numericalDerivative11<Vector3, Pose3>(
      boost::bind(rkfErr, rkf1, PWtoBi, _1, PWtoBj, PGjtoBj, tGtoB_G),
        PGitoBi, delta);
  Matrix H3e = numericalDerivative11<Vector3, Pose3>(
      boost::bind(rkfErr, rkf1, PWtoBi, PGitoBi, _1, PGjtoBj, tGtoB_G),
        PWtoBj, delta);
  Matrix H4e = numericalDerivative11<Vector3, Pose3>(
      boost::bind(rkfErr, rkf1, PWtoBi, PGitoBi, PWtoBj, _1, tGtoB_G),
        PGjtoBj, delta);
  Matrix H5e = numericalDerivative11<Vector3, Point3>(
      boost::bind(rkfErr, rkf1, PWtoBi, PGitoBi, PWtoBj, PGjtoBj, _1),
      tGtoB_G, delta);
  EXPECT(assert_equal(H1e,H1,1e-10));
  EXPECT(assert_equal(H2e,H2,1e-10));
  EXPECT(assert_equal(H3e,H3,1e-10));
  EXPECT(assert_equal(H4e,H4,1e-10));
  EXPECT(assert_equal(H5e,H5,1e-10));

}

/* ************************************************************************* */
int main () {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */


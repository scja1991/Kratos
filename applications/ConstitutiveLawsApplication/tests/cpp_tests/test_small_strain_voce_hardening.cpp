// License: BSD License (Kratos default license: kratos/license.txt)
#include "custom_constitutive/auxiliary_files/hardening_laws/voce_hardening_law.h"
#include <iomanip>
#include <sstream>
#include "containers/model.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/table_accessor.h"
#include "tests/cpp_tests/constitutive_laws_fast_suite.h"
#include "tests/cpp_tests/j2_hardening_test_utilities.h"

namespace Kratos::Testing
{
namespace VoceTests
{
using namespace J2HardeningTestUtilities;
class Inspectable : public SmallStrainJ2IsotropicHardening3D
{
public:
    using SmallStrainJ2IsotropicHardening3D::EvaluateTrial;
    double Cache() { return GetThreshold(); }
};
struct Point
{
    Model ModelInstance;
    ModelPart& Part;
    Tetrahedra3D4<Node>::Pointer Geometry;
    Properties properties;
    Vector strain = ZeroVector(6), stress = ZeroVector(6), N = Vector(4,0.25);
    Matrix tangent = ZeroMatrix(6,6);
    ConstitutiveLaw::Parameters values;
    Inspectable law;
    Point() : Part(ModelInstance.CreateModelPart("Voce"))
    {
        Geometry = Kratos::make_shared<Tetrahedra3D4<Node>>(
            Part.CreateNewNode(1,0,0,0),Part.CreateNewNode(2,1,0,0),
            Part.CreateNewNode(3,0,1,0),Part.CreateNewNode(4,0,0,1));
        properties[YOUNG_MODULUS]=210000.0; properties[POISSON_RATIO]=0.3;
        properties[HARDENING_CURVE]=8; properties[YIELD_STRESS]=250.0;
        properties[EXPONENTIAL_SATURATION_YIELD_STRESS]=450.0; properties[HARDENING_EXPONENT]=20.0;
        values.SetElementGeometry(*Geometry); values.SetMaterialProperties(properties);
        values.SetProcessInfo(Part.GetProcessInfo()); values.SetShapeFunctionsValues(N);
        values.SetStrainVector(strain); values.SetStressVector(stress); values.SetConstitutiveMatrix(tangent);
        values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        Reset();
    }
    void Reset() { law.ResetMaterial(properties,*Geometry,N); }
    void Evaluate() { law.CalculateMaterialResponseCauchy(values); }
    void Commit() { law.FinalizeMaterialResponseCauchy(values); }
    double P(ConstitutiveLaw* pLaw=nullptr)
    { double value=0.0; (pLaw?pLaw:&law)->GetValue(ACCUMULATED_PLASTIC_STRAIN,value); return value; }
    Vector Plastic(ConstitutiveLaw* pLaw=nullptr)
    { Vector value; (pLaw?pLaw:&law)->GetValue(PLASTIC_STRAIN_VECTOR,value); return value; }
};
// Independent long-double formula, including products outside double range.
long double yield_long(const Properties& rProperties, long double P)
{
    const long double initial=rProperties[YIELD_STRESS];
    const long double q=static_cast<long double>(rProperties[EXPONENTIAL_SATURATION_YIELD_STRESS])-initial;
    const long double b=rProperties[HARDENING_EXPONENT];
    return initial-q*std::expm1(-b*P);
}
double yield_reference(const Properties& rProperties, double P)
{ return static_cast<double>(yield_long(rProperties,P)); }
struct Reference { Vector Stress, Plastic; double P; };
Reference reference_step(Point& rPoint)
{
    const Matrix c=elastic_reference(rPoint.properties);
    const Vector old=rPoint.Plastic(), predictor=prod(c,rPoint.strain-old);
    const double q=q_reference(predictor), old_p=rPoint.P();
    if (q<=yield_reference(rPoint.properties,old_p)) return {predictor,old,old_p};
    const double three_g=1.5*rPoint.properties[YOUNG_MODULUS]/(1.0+rPoint.properties[POISSON_RATIO]);
    double lower=0.0,upper=(q-yield_reference(rPoint.properties,old_p))/three_g;
    for(int i=0;i<120;++i) {
        const double mid=(lower+upper)/2.0;
        if(q-three_g*mid-yield_reference(rPoint.properties,old_p+mid)>0.0) lower=mid;
        else upper=mid;
    }
    const double lambda=(lower+upper)/2.0;
    Matrix flow_tensor=MathUtils<double>::StressVectorToTensor(predictor);
    const double mean=(flow_tensor(0,0)+flow_tensor(1,1)+flow_tensor(2,2))/3.0;
    for(SizeType i=0;i<3;++i) flow_tensor(i,i)-=mean;
    flow_tensor*=1.5/q;
    const Vector increment=lambda*MathUtils<double>::StrainTensorToVector(flow_tensor,6);
    return {Vector(predictor-prod(c,increment)),Vector(old+increment),old_p+lambda};
}
void compare(Point& rPoint,const Reference& rReference)
{
    rPoint.Evaluate();
    const auto result=rPoint.law.EvaluateTrial(rPoint.values,true);
    KRATOS_EXPECT_TRUE(result.Candidate.has_value());
    KRATOS_EXPECT_LE(result.StressResidual,result.StressTolerance);
    KRATOS_EXPECT_LE(result.YieldResidual,result.StressTolerance);
    if(result.Status==SmallStrainJ2IsotropicHardening3D::Integrator::IntegrationStatus::Converged) {
        KRATOS_EXPECT_LE(std::abs(result.YieldResidual),result.StressTolerance);
    }
    KRATOS_EXPECT_VECTOR_NEAR(rPoint.stress,rReference.Stress,1.0e-7+1.0e-9*norm_inf(rReference.Stress));
    KRATOS_EXPECT_NEAR(result.Candidate->P,rReference.P,1.0e-12+1.0e-9*rReference.P);
    KRATOS_EXPECT_VECTOR_NEAR(result.Candidate->PlasticStrain,rReference.Plastic,1.0e-12);
}
} // namespace VoceTests

KRATOS_TEST_CASE_IN_SUITE(VoceEndpointScalarFunction,KratosConstitutiveLawsFastSuite)
{
    VoceTests::Point point;
    double previous=0.0;
    for(double p:{0.0,1.0e-12,1.0e-8,0.002,0.01,0.1,1.0,10.0}) {
        double y=0.0,h=0.0;VoceHardeningLaw::CalculateStressThreshold(point.properties,p,y,h);
        const long double ref=VoceTests::yield_long(point.properties,p);
        const long double href=200.0L*20.0L*std::exp(-20.0L*p);
        KRATOS_EXPECT_NEAR(y,static_cast<double>(ref),1.0e-12*y);
        KRATOS_EXPECT_NEAR(h,static_cast<double>(href),1.0e-12*std::max(1.0,h));
        KRATOS_EXPECT_GE(y,previous);KRATOS_EXPECT_LE(y,450.0);previous=y;
        if(p>1.0e-6 && p<1.0) {
            const long double step=5.0e-8L;
            const double difference=static_cast<double>((VoceTests::yield_long(point.properties,p+step)-
                VoceTests::yield_long(point.properties,p-step))/(2.0L*step));
            KRATOS_EXPECT_NEAR(h,difference,1.0e-8*std::max(1.0,h));
        }
    }
    double y=0.0,h=0.0;VoceHardeningLaw::CalculateStressThreshold(point.properties,0.0,y,h);
    KRATOS_EXPECT_EQ(y,250.0);KRATOS_EXPECT_EQ(h,4000.0);
    VoceHardeningLaw::CalculateStressThreshold(point.properties,1000.0,y,h);
    KRATOS_EXPECT_EQ(y,450.0);KRATOS_EXPECT_EQ(h,0.0);
    for (bool zero_q:{false,true}) {
        VoceTests::Point limit;
        if(zero_q) limit.properties[EXPONENTIAL_SATURATION_YIELD_STRESS]=250.0;
        else limit.properties[HARDENING_EXPONENT]=0.0;
        KRATOS_EXPECT_EQ(VoceHardeningLaw::Check(limit.properties),0);
        for(double p:{0.0,0.01,1.0e300}) {
            VoceHardeningLaw::CalculateStressThreshold(limit.properties,p,y,h);
            KRATOS_EXPECT_EQ(y,250.0);KRATOS_EXPECT_EQ(h,0.0);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(VoceEndpointExtremeArguments,KratosConstitutiveLawsFastSuite)
{
    VoceTests::Point point;auto& props=point.properties;
    props[YIELD_STRESS]=1.0;props[EXPONENTIAL_SATURATION_YIELD_STRESS]=1.0e300;props[HARDENING_EXPONENT]=1.0e-300;
    double y=0.0,h=0.0;
    VoceHardeningLaw::CalculateStressThreshold(props,1.0e-10,y,h);
    KRATOS_EXPECT_NEAR(y,1.0+1.0e-10,1.0e-14);KRATOS_EXPECT_NEAR(h,1.0,1.0e-14);
    props[YIELD_STRESS]=1.0e-308;props[EXPONENTIAL_SATURATION_YIELD_STRESS]=1.0e308;
    for(double p:{1.0e-20,1.0e-300}) {
        VoceHardeningLaw::CalculateStressThreshold(props,p,y,h);
        const double ref=1.0e-308+1.0e8*p; // 1-exp(-bp)=bp to double resolution here.
        KRATOS_EXPECT_GT(y,1.0e-308);KRATOS_EXPECT_NEAR(y,ref,1.0e-12*ref);
    }
    props[YIELD_STRESS]=1.0;
    props[EXPONENTIAL_SATURATION_YIELD_STRESS]=1.0e308;props[HARDENING_EXPONENT]=1.0;
    VoceHardeningLaw::CalculateStressThreshold(props,800.0,y,h);
    // 80-digit Decimal references, rounded here; no extended long-double range required.
    const double reference=3.667874584177687e-40;
    KRATOS_EXPECT_GT(h,0.0);KRATOS_EXPECT_NEAR(h,reference,1.0e-12*reference);
    props[HARDENING_EXPONENT]=10.0;
    for(const auto& sample:{std::make_pair(1.0,4.539992976248485e+304),
                            std::make_pair(100.0,5.075958897549457e-126)}) {
        KRATOS_EXPECT_TRUE(VoceHardeningLaw::TryCalculateStressThreshold(props,sample.first,y,h));
        KRATOS_EXPECT_NEAR(h,sample.second,1.0e-12*sample.second);
    }
    KRATOS_EXPECT_FALSE(VoceHardeningLaw::TryCalculateStressThreshold(props,0.0,y,h));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(VoceHardeningLaw::Check(props),"not representable");
    props[EXPONENTIAL_SATURATION_YIELD_STRESS]=2.0;
    props[HARDENING_EXPONENT]=std::numeric_limits<double>::max();
    VoceHardeningLaw::CalculateStressThreshold(props,std::numeric_limits<double>::max(),y,h);
    KRATOS_EXPECT_EQ(y,2.0);KRATOS_EXPECT_EQ(h,0.0);
}

KRATOS_TEST_CASE_IN_SUITE(VoceEndpointParametersAndSelection,KratosConstitutiveLawsFastSuite)
{
    VoceTests::Point point;const auto original=point.properties;
    KRATOS_EXPECT_FALSE(point.properties.Has(SWIFT_COEFFICIENT));
    KRATOS_EXPECT_FALSE(point.properties.Has(FRACTURE_ENERGY));
    const std::string before=J2HardeningTestUtilities::state_trace(point.law);
    auto check=[&](){return point.law.Check(point.properties,*point.Geometry,point.Part.GetProcessInfo());};
    for(const auto* variable:{&YIELD_STRESS,&EXPONENTIAL_SATURATION_YIELD_STRESS,&HARDENING_EXPONENT}) {
        point.properties=original;point.properties.GetData().Erase(*variable);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"Missing");
        for(double value:{-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
            point.properties=original;point.properties[*variable]=value;
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"");
            double y=0.0,h=0.0;KRATOS_EXPECT_FALSE(VoceHardeningLaw::TryCalculateStressThreshold(point.properties,0.01,y,h));
        }
    }
    point.properties=original;point.properties[YIELD_STRESS]=0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"Invalid Voce");
    point.properties=original;point.properties[EXPONENTIAL_SATURATION_YIELD_STRESS]=249.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"Invalid Voce");
    for(double value:{1.0,-1.0,std::numeric_limits<double>::quiet_NaN()}) {
        point.properties=original;point.properties[ISOTROPIC_HARDENING_MODULUS]=value;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"ISOTROPIC_HARDENING_MODULUS");
    }
    point.properties=original;point.properties[ISOTROPIC_HARDENING_MODULUS]=0.0;KRATOS_EXPECT_EQ(check(),0);
    for(const auto* variable:{&YIELD_STRESS_TENSION,&YIELD_STRESS_COMPRESSION}) {
        point.properties=original;point.properties[*variable]=250.0;KRATOS_EXPECT_EQ(check(),0);
        point.properties[*variable]=251.0;KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"Redundant yield stress");
    }
    for(const auto* variable:{&YIELD_STRESS,&EXPONENTIAL_SATURATION_YIELD_STRESS,&HARDENING_EXPONENT}) {
        // Fresh fixtures: Properties assignment does not remove existing accessors.
        VoceTests::Point table_point;
        Table<double,double> table;table.PushBack(0.0,1.0);table.PushBack(100.0,2.0);
        table_point.properties.SetTable(TEMPERATURE,*variable,table);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(table_point.law.Check(table_point.properties,*table_point.Geometry,
            table_point.Part.GetProcessInfo()),"constant isothermal");
        VoceTests::Point accessor_point;
        TableAccessor accessor(TEMPERATURE,"node_historical");accessor_point.properties.SetAccessor(*variable,accessor.Clone());
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(accessor_point.law.Check(accessor_point.properties,*accessor_point.Geometry,
            accessor_point.Part.GetProcessInfo()),"constant isothermal");
    }
    point.properties=original;point.properties[THERMAL_EXPANSION_COEFFICIENT]=1.0e-5;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"Thermal J2");
    point.properties=original;point.properties[HARDENING_CURVE]=9;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"Unsupported HARDENING_CURVE");
    point.properties=original;point.properties.GetData().Erase(HARDENING_CURVE);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(),"requires HARDENING_CURVE");
    point.properties=original;SmallStrainSwiftJ2Plasticity3D swift;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(swift.Check(point.properties,*point.Geometry,point.Part.GetProcessInfo()),"conflicting HARDENING_CURVE");
    point.properties=original;
    for(double p:{-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        double y=0.0,h=0.0;KRATOS_EXPECT_FALSE(VoceHardeningLaw::TryCalculateStressThreshold(point.properties,p,y,h));
    }
    KRATOS_EXPECT_EQ(J2HardeningTestUtilities::state_trace(point.law),before);
}

KRATOS_TEST_CASE_IN_SUITE(VoceEndpointPathsAndHistory,KratosConstitutiveLawsFastSuite)
{
    using namespace J2HardeningTestUtilities;
    for(int mode:{0,1,2}) {
        VoceTests::Point point;
        if(mode==1) point.properties[EXPONENTIAL_SATURATION_YIELD_STRESS]=250.0;
        if(mode==2) point.properties[HARDENING_EXPONENT]=0.0;
        point.Reset();
        for(const Vector& strain:{make_vector({0.0001,0,0,0,0,0}),make_vector({0.006,0,0,0,0,0}),
            make_vector({0.005,0,0,0,0,0}),make_vector({0.008,0,0,0,0,0}),make_vector({-0.008,0,0,0,0,0}),
            make_vector({0,0,0,0.02,0,0}),make_vector({0.006,-0.004,-0.002,-0.012,0.009,0.004})}) {
            point.strain=strain;const auto reference=VoceTests::reference_step(point);
            const std::string state=state_trace(point.law);
            VoceTests::compare(point,reference);
            KRATOS_EXPECT_EQ(state_trace(point.law),state);
            point.Commit();KRATOS_EXPECT_NEAR(point.P(),reference.P,1.0e-12);
            KRATOS_EXPECT_VECTOR_NEAR(point.Plastic(),reference.Plastic,1.0e-12);
            const double p=point.P();point.Commit();KRATOS_EXPECT_EQ(point.P(),p);
        }
        KRATOS_EXPECT_GT(point.P(),AdvancedConstitutiveLawUtilities<6>::CalculateEquivalentPlasticStrainIncrement(point.Plastic())+0.001);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VoceEndpointLifecycleAndRestart,KratosConstitutiveLawsFastSuite)
{
    using namespace J2HardeningTestUtilities;VoceTests::Point point;
    point.strain=make_vector({0.007,-0.001,0,0.004,0.002,0});point.Commit();
    VoceTests::Inspectable copy(point.law),assigned;assigned=point.law;
    auto clone=point.law.Clone();KRATOS_EXPECT_EQ(clone->Info(),"SmallStrainJ2IsotropicHardening3DLaw");
    KRATOS_EXPECT_TRUE(typeid(*clone)==typeid(SmallStrainJ2IsotropicHardening3D));
    VoceTests::Inspectable imported;imported.SetPlasticState(point.Plastic(),point.P());
    KRATOS_EXPECT_EQ(imported.Cache(),0.0);
    const std::string imported_state=state_trace(imported);double yield=0.0;
    imported.CalculateValue(point.values,YIELD_STRESS,yield);
    KRATOS_EXPECT_NEAR(yield,VoceTests::yield_reference(point.properties,point.P()),1.0e-12*yield);
    KRATOS_EXPECT_EQ(state_trace(imported),imported_state);
    for(auto* law:{static_cast<ConstitutiveLaw*>(&copy),static_cast<ConstitutiveLaw*>(&assigned),clone.get()}) {
        KRATOS_EXPECT_EQ(point.P(law),point.P());KRATOS_EXPECT_VECTOR_NEAR(point.Plastic(law),point.Plastic(),0.0);
    }
    for(auto trace:{Serializer::SERIALIZER_NO_TRACE,Serializer::SERIALIZER_TRACE_ERROR}) {
        StreamSerializer serializer(trace);serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
        serializer.save("Law",point.law);serializer.save("Sentinel",42.125);
        VoceTests::Inspectable restored;serializer.load("Law",restored);
        double sentinel=0.0;serializer.load("Sentinel",sentinel);KRATOS_EXPECT_EQ(sentinel,42.125);
        KRATOS_EXPECT_EQ(restored.Cache(),0.0);
        const std::string state=state_trace(restored);restored.CalculateValue(point.values,YIELD_STRESS,yield);
        KRATOS_EXPECT_EQ(state_trace(restored),state);
        KRATOS_EXPECT_NEAR(yield,VoceTests::yield_reference(point.properties,point.P()),1.0e-12*yield);
        auto original=point.law.Clone();
        for(const Vector& strain:{make_vector({0.004,0,0,0.002,0,0}),make_vector({0.012,-0.002,0,-0.007,0.005,0})}) {
            point.strain=strain;original->CalculateMaterialResponseCauchy(point.values);
            const Vector stress=point.stress;const Matrix tangent=point.tangent;
            restored.CalculateMaterialResponseCauchy(point.values);
            KRATOS_EXPECT_VECTOR_NEAR(point.stress,stress,1.0e-7);KRATOS_EXPECT_MATRIX_NEAR(point.tangent,tangent,1.0e-5);
            original->FinalizeMaterialResponseCauchy(point.values);restored.FinalizeMaterialResponseCauchy(point.values);
            KRATOS_EXPECT_NEAR(point.P(original.get()),point.P(&restored),1.0e-14);
        }
    }
    point.strain=make_vector({0.03,-0.005,0.002,0.025,0.003,0});
    point.properties[MAX_NUMBER_NL_CL_ITERATIONS]=1;
    point.stress=ScalarVector(6,-123.0);point.tangent=ScalarMatrix(6,6,-456.0);
    const std::string state=state_trace(point.law);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"iteration limit");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Commit(),"iteration limit");
    KRATOS_EXPECT_EQ(state_trace(point.law),state);
    KRATOS_EXPECT_VECTOR_NEAR(point.stress,Vector(ScalarVector(6,-123.0)),0.0);
    KRATOS_EXPECT_MATRIX_NEAR(point.tangent,Matrix(ScalarMatrix(6,6,-456.0)),0.0);
    point.properties[MAX_NUMBER_NL_CL_ITERATIONS]=100;
    VoceTests::compare(point,VoceTests::reference_step(point));point.Commit();
    point.Reset();KRATOS_EXPECT_EQ(point.P(),0.0);KRATOS_EXPECT_EQ(norm_2(point.Plastic()),0.0);KRATOS_EXPECT_EQ(point.law.Cache(),250.0);
}

KRATOS_TEST_CASE_IN_SUITE(VoceEndpointDiscreteTangent,KratosConstitutiveLawsFastSuite)
{
    using namespace J2HardeningTestUtilities;
    const std::array<const char*,9> names{{"elastic","tension","compression","shear","unload","reload","nonproportional","below_yield","above_yield"}};
    double maximum_matrix=0.0,maximum_direction=0.0;int count=0,accepted=0,faithful=0,minimum_adjacent=9;
    auto record=[&](VoceTests::Point& point,bool plastic,const std::string& name) {
        SCOPED_TRACE(name);
        VoceTests::compare(point,VoceTests::reference_step(point));
        const auto result=check_tangent_sweep(point,plastic,VoceTests::yield_reference);
        ++count;accepted+=result.AcceptedSteps;faithful+=result.BranchFaithfulSteps;
        minimum_adjacent=std::min(minimum_adjacent,result.MaxConsecutive);
        maximum_matrix=std::max(maximum_matrix,result.MaxMatrixError);maximum_direction=std::max(maximum_direction,result.MaxDirectionError);
        std::ostringstream data;data<<std::setprecision(17)<<"branch_faithful="<<result.BranchFaithfulSteps
            <<";accepted="<<result.AcceptedSteps<<";adjacent="<<result.MaxConsecutive<<";max_matrix="<<result.MaxMatrixError
            <<";max_direction="<<result.MaxDirectionError;::testing::Test::RecordProperty(name,data.str());
    };
    for(int mode:{0,1,2}) for(int path=0;path<9;++path) {
        VoceTests::Point point;
        if(mode==1) point.properties[EXPONENTIAL_SATURATION_YIELD_STRESS]=250.0;
        if(mode==2) point.properties[HARDENING_EXPONENT]=0.0;
        point.Reset();bool plastic=true;const double e=210000.0,nu=0.3,two_g=e/(1.0+nu);
        switch(path) {
            case 0:point.strain[0]=0.0001;plastic=false;break;
            case 1:case 2:{
                double lower=0.0,upper=0.008;
                for(int i=0;i<120;++i) {const double mid=(lower+upper)/2.0;
                    if(e*(0.008-mid)>VoceTests::yield_reference(point.properties,mid)) lower=mid;else upper=mid;}
                const double p=(lower+upper)/2.0,sigma=e*(0.008-p),sign=path==1?1.0:-1.0;
                point.strain=make_vector({sign*0.008,sign*(-nu*sigma/e-p/2.0),sign*(-nu*sigma/e-p/2.0),0,0,0});
                point.Evaluate();
                KRATOS_EXPECT_VECTOR_NEAR(point.stress,make_vector({sign*sigma,0,0,0,0,0}),1.0e-7);
                const auto trial=point.law.EvaluateTrial(point.values,false);
                KRATOS_EXPECT_NEAR(trial.Candidate->P,p,1.0e-12);break;
            }
            case 3:point.strain[3]=0.01;break;
            case 4:case 5:
                point.strain[0]=0.008;point.Commit();point.strain=point.Plastic();
                point.strain[0]+=0.5*VoceTests::yield_reference(point.properties,point.P())/two_g;
                if(path==4) plastic=false;else {point.Commit();point.strain=make_vector({0.012,0,0,0,0,0});}break;
            case 6:point.strain=make_vector({0.007,-0.001,-0.001,0.003,0,0});point.Commit();
                point.strain=make_vector({0.012,-0.004,0.001,-0.006,0.004,0.002});break;
            case 7:case 8:plastic=path==8;point.strain[0]=250.0/two_g*(1.0+(plastic?1.0:-1.0)*1.0e-4);break;
        }
        record(point,plastic,"mode"+std::to_string(mode)+"_"+names[path]);
    }
    VoceTests::Point saturation;saturation.law.SetPlasticState(ZeroVector(6),0.5);saturation.strain[0]=0.008;
    record(saturation,true,"saturation_bp10");
    ::testing::Test::RecordProperty("tangent_cases",count);::testing::Test::RecordProperty("accepted_steps",accepted);
    ::testing::Test::RecordProperty("branch_faithful_steps",faithful);::testing::Test::RecordProperty("minimum_adjacent_steps",minimum_adjacent);
    for(const auto& metric:{std::make_pair("max_relative_matrix_error",maximum_matrix),std::make_pair("max_relative_direction_error",maximum_direction)}) {
        std::ostringstream value;value<<std::setprecision(17)<<metric.second;::testing::Test::RecordProperty(metric.first,value.str());
    }
}

KRATOS_TEST_CASE_IN_SUITE(VoceEndpointRepresentabilityTransitions,KratosConstitutiveLawsFastSuite)
{
    // 400-digit Decimal.from_float references: every input is an exact binary
    // double. q = saturation - initial is formed at that precision before exp.
    // Hex literals freeze those inputs and the independently rounded references.
    struct Case { const char* Name; double Initial, Saturation, B, P; bool Finite; double Yield, H; };
    const Case cases[] = {
        {"review", 0x1.0000000000000p+0, 0x1.1ccf385ebc8a0p+1023, 0x1.4000000000000p+3, 0x1.5f740f2a2d267p-3, false, 0x1.d3380a5712aa5p+1022, 0.0},
        {"near_max_overflow", 0x1.0000000000000p+0, 0x1.1ccf385ebc8a0p+1023, 0x1.4000000000000p+3, 0x1.5f740f2a2d2ccp-3, false, 0x1.d3380a5712ad8p+1022, 0.0},
        {"near_max_finite", 0x1.0000000000000p+0, 0x1.1ccf385ebc8a0p+1023, 0x1.4000000000000p+3, 0x1.5f740f2a2d2d4p-3, true, 0x1.d3380a5712adcp+1022, 0x1.ffffffffffff4p+1023},
        {"min_below", 0x0.730d67819e8d2p-1022, 0x1.1ccf385ebc8a0p+1023, 0x1.0000000000000p+0, 0x0.fffffffffffffp-1022, true, 0x1.1ccf385ebc89fp+1, 0x1.1ccf385ebc8a0p+1023},
        {"min_above", 0x0.730d67819e8d2p-1022, 0x1.1ccf385ebc8a0p+1023, 0x1.0000000000000p+0, 0x1.0000000000001p-1022, true, 0x1.1ccf385ebc8a1p+1, 0x1.1ccf385ebc8a0p+1023},
        {"log2_below", 0x1.f400000000000p+7, 0x1.c200000000000p+8, 0x1.0000000000000p+0, 0x1.62e42fefa39eep-1, true, 0x1.5e00000000000p+8, 0x1.9000000000001p+6},
        {"log2_above", 0x1.f400000000000p+7, 0x1.c200000000000p+8, 0x1.0000000000000p+0, 0x1.62e42fefa39f0p-1, true, 0x1.5e00000000000p+8, 0x1.8ffffffffffffp+6},
        {"exp_normal_below", 0x1.0000000000000p+0, 0x1.1ccf385ebc8a0p+1023, 0x1.0000000000000p+0, 0x1.6232bdd7abcd1p+9, true, 0x1.1ccf385ebc8a0p+1023, 0x1.1ccf385ebcb63p+1},
        {"exp_normal_above", 0x1.0000000000000p+0, 0x1.1ccf385ebc8a0p+1023, 0x1.0000000000000p+0, 0x1.6232bdd7abcd3p+9, true, 0x1.1ccf385ebc8a0p+1023, 0x1.1ccf385ebc6f0p+1},
        {"decay_140.0", 0x1.0000000000000p+0, 0x1.1ccf385ebc8a0p+1023, 0x1.4000000000000p+3, 0x1.1800000000000p+7, true, 0x1.1ccf385ebc8a0p+1023, 0x1.a0a8c247c1a32p-994},
        {"decay_145.0", 0x1.0000000000000p+0, 0x1.1ccf385ebc8a0p+1023, 0x1.4000000000000p+3, 0x1.2200000000000p+7, true, 0x1.1ccf385ebc8a0p+1023, 0x0.000000000017cp-1022},
        {"decay_150.0", 0x1.0000000000000p+0, 0x1.1ccf385ebc8a0p+1023, 0x1.4000000000000p+3, 0x1.2c00000000000p+7, true, 0x1.1ccf385ebc8a0p+1023, 0x0.0p+0},
    };
    VoceTests::Point point;
    double maximum_yield_error=0.0,maximum_h_error=0.0;
    for(const auto& sample:cases) {
        SCOPED_TRACE(sample.Name);
        point.properties[YIELD_STRESS]=sample.Initial;
        point.properties[EXPONENTIAL_SATURATION_YIELD_STRESS]=sample.Saturation;
        point.properties[HARDENING_EXPONENT]=sample.B;
        double y=-123.0,h=-456.0;
        const bool success=VoceHardeningLaw::TryCalculateStressThreshold(point.properties,sample.P,y,h);
        KRATOS_EXPECT_EQ(success,sample.Finite);
        if(!sample.Finite) {
            KRATOS_EXPECT_EQ(y,-123.0);KRATOS_EXPECT_EQ(h,-456.0);
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(VoceHardeningLaw::CalculateStressThreshold(
                point.properties,sample.P,y,h),"not representable");
            KRATOS_EXPECT_EQ(y,-123.0);KRATOS_EXPECT_EQ(h,-456.0);
        } else {
            // 8/32 eps allow ordinary exp and a few scaled products, without
            // hiding the status decision. Two subnormal quanta cover rounding.
            const double yield_tolerance=std::max(8.0*std::numeric_limits<double>::epsilon()*sample.Yield,
                2.0*std::numeric_limits<double>::denorm_min());
            const double h_tolerance=std::max(32.0*std::numeric_limits<double>::epsilon()*sample.H,
                2.0*std::numeric_limits<double>::denorm_min());
            KRATOS_EXPECT_NEAR(y,sample.Yield,yield_tolerance);
            KRATOS_EXPECT_NEAR(h,sample.H,h_tolerance);
            maximum_yield_error=std::max(maximum_yield_error,std::abs(y-sample.Yield)/sample.Yield);
            if(sample.H>=std::numeric_limits<double>::min())
                maximum_h_error=std::max(maximum_h_error,std::abs(h-sample.H)/sample.H);
        }
    }
    // The same parameters can admit an individual p>0 evaluation while the
    // material is inadmissible because its initial derivative overflows.
    point.properties[YIELD_STRESS]=1.0;
    point.properties[EXPONENTIAL_SATURATION_YIELD_STRESS]=1.0e308;
    point.properties[HARDENING_EXPONENT]=10.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(VoceHardeningLaw::Check(point.properties),"not representable");
    double y=-123.0,h=-456.0;
    KRATOS_EXPECT_FALSE(VoceHardeningLaw::TryCalculateStressThreshold(point.properties,0.0,y,h));
    KRATOS_EXPECT_EQ(y,-123.0);KRATOS_EXPECT_EQ(h,-456.0);
    for(const auto& metric:{std::make_pair("max_yield_relative_error",maximum_yield_error),
        std::make_pair("max_normal_h_relative_error",maximum_h_error)}) {
        std::ostringstream value;value<<std::setprecision(17)<<metric.second;
        ::testing::Test::RecordProperty(metric.first,value.str());
    }
}

KRATOS_TEST_CASE_IN_SUITE(VoceEndpointParameterIsolation,KratosConstitutiveLawsFastSuite)
{
    auto prepare_swift=[](VoceTests::Point& point) {
        point.properties[HARDENING_CURVE]=7;
        point.properties.GetData().Erase(YIELD_STRESS);
        point.properties[SWIFT_COEFFICIENT]=1000.0;
        point.properties[SWIFT_INITIAL_STRAIN]=0.01;
        point.properties[SWIFT_HARDENING_EXPONENT]=0.2;
    };
    auto check=[](VoceTests::Point& point) {
        return point.law.Check(point.properties,*point.Geometry,point.Part.GetProcessInfo());
    };
    Table<double,double> table;table.PushBack(0.0,1.0);table.PushBack(100.0,2.0);
    TableAccessor accessor(TEMPERATURE,"node_historical");
    for(bool swift:{false,true}) {
        VoceTests::Point point;
        if(swift) prepare_swift(point);
        point.Reset();point.strain[0]=0.006;point.Evaluate();
        const Vector stress=point.stress;const Matrix tangent=point.tangent;
        if(swift) {
            // Initial YIELD_STRESS is redundant in Swift, not an inactive field.
            // Only saturation and b are pure inactive Voce fields here.
            point.properties[EXPONENTIAL_SATURATION_YIELD_STRESS]=-1.0;
            point.properties[HARDENING_EXPONENT]=std::numeric_limits<double>::quiet_NaN();
            point.properties.SetTable(TEMPERATURE,EXPONENTIAL_SATURATION_YIELD_STRESS,table);
            point.properties.SetAccessor(HARDENING_EXPONENT,accessor.Clone());
            SmallStrainSwiftJ2Plasticity3D wrapper;
            KRATOS_EXPECT_EQ(wrapper.Check(point.properties,*point.Geometry,point.Part.GetProcessInfo()),0);
        } else {
            point.properties[SWIFT_COEFFICIENT]=std::numeric_limits<double>::quiet_NaN();
            point.properties[SWIFT_INITIAL_STRAIN]=-1.0;
            point.properties[SWIFT_HARDENING_EXPONENT]=std::numeric_limits<double>::infinity();
            point.properties.SetTable(TEMPERATURE,SWIFT_COEFFICIENT,table);
            point.properties.SetAccessor(SWIFT_INITIAL_STRAIN,accessor.Clone());
        }
        KRATOS_EXPECT_EQ(check(point),0);point.Evaluate();
        KRATOS_EXPECT_VECTOR_NEAR(point.stress,stress,0.0);
        KRATOS_EXPECT_MATRIX_NEAR(point.tangent,tangent,0.0);
    }
    for(bool swift:{false,true}) {
        for(const auto* variable:{&ISOTROPIC_HARDENING_MODULUS,&YIELD_STRESS,&YIELD_STRESS_TENSION,&YIELD_STRESS_COMPRESSION}) {
            if(swift && variable==&ISOTROPIC_HARDENING_MODULUS) continue; // No Swift contract for this inactive field.
            for(bool use_table:{false,true}) {
                VoceTests::Point point;if(swift) prepare_swift(point);
                const double initial=swift?1000.0*std::pow(0.01,0.2):250.0;
                point.properties[*variable]=variable==&ISOTROPIC_HARDENING_MODULUS?0.0:initial;
                KRATOS_EXPECT_EQ(check(point),0);
                if(use_table) point.properties.SetTable(TEMPERATURE,*variable,table);
                else point.properties.SetAccessor(*variable,accessor.Clone());
                KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(point),"constant isothermal");
            }
        }
        // Redundant yields in Swift remain checked, even though saturation/b
        // are inactive. They cannot be reclassified as irrelevant Voce fields.
        for(const auto* variable:{&YIELD_STRESS,&YIELD_STRESS_TENSION,&YIELD_STRESS_COMPRESSION}) {
            if(!swift && variable==&YIELD_STRESS) continue; // Voce's primary parameter, not redundant.
            VoceTests::Point point;if(swift) prepare_swift(point);
            const double initial=swift?1000.0*std::pow(0.01,0.2):250.0;
            point.properties[*variable]=initial+1.0;
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(check(point),"Redundant yield stress");
        }
    }
}
} // namespace Kratos::Testing

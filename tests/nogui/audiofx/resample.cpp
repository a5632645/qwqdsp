#include "AudioFile.h"
#include "work_dir.hpp"
#include <qwqdsp/filter/window_fir.hpp>
#include <qwqdsp/fx/polyphase_resample_fir.hpp>
#include <qwqdsp/fx/resample.hpp>
#include <qwqdsp/fx/resample_coeffs.h>
#include <qwqdsp/fx/resample_iir.hpp>
#include <qwqdsp/window/blackman_nuttall.hpp>


constexpr float kTargetFs = 88200.0f;

static void ResampleIir() {
    AudioFile<float> infile;
    infile.load(qwqdsp_support::InputFile("sweep.wav"));
    auto& sweep = infile.samples.front();

    qwqdsp_fx::ResampleIIR<qwqdsp_fx::coeff::BestCoeffs<float>, 128> resample;
    resample.Init(infile.getSampleRate(), kTargetFs);
    auto sweep_resample = resample.Process<float>(sweep);

    AudioFile<float> outfile;
    outfile.setNumChannels(1);
    outfile.setBitDepth(32);
    outfile.setSampleRate(kTargetFs);
    outfile.setNumSamplesPerChannel(sweep_resample.size());
    outfile.samples.front() = std::move(sweep_resample);
    outfile.save(qwqdsp_support::OutputFile("sweep_resample_iir.wav"));
}

static void ResampleFir() {
    AudioFile<float> infile;
    infile.load(qwqdsp_support::InputFile("sweep.wav"));
    auto& sweep = infile.samples.front();

    qwqdsp_fx::Resample resample;
    resample.Init(infile.getSampleRate(), kTargetFs, 180, 127, 128);
    auto sweep_resample = resample.Process(sweep);

    AudioFile<float> outfile;
    outfile.setNumChannels(1);
    outfile.setBitDepth(32);
    outfile.setSampleRate(kTargetFs);
    outfile.setNumSamplesPerChannel(sweep_resample.size());
    outfile.samples.front() = std::move(sweep_resample);
    outfile.save(qwqdsp_support::OutputFile("sweep_resample_fir.wav"));
}

static void ResamplePolyphaseFir() {
    constexpr int coeff_len = 1024;
    constexpr int ratio = 2;
    float coeff[coeff_len];
    qwqdsp_filter::WindowFIR::Lowpass(coeff, std::numbers::pi_v<float> / ratio);
    qwqdsp_window::BlackmanNuttall::ApplyWindow(coeff, false);

    AudioFile<float> infile;
    infile.load(qwqdsp_support::InputFile("sweep.wav"));
    auto& sweep = infile.samples.front();

    {
        qwqdsp_fx::PolyphaseDownsamplerFir resample;
        resample.Init(coeff, ratio);
        resample.Reset();
        auto sweep_resample = resample.Process(sweep);

        AudioFile<float> outfile;
        outfile.setNumChannels(1);
        outfile.setBitDepth(32);
        outfile.setSampleRate(infile.getSampleRate() / ratio);
        outfile.setNumSamplesPerChannel(sweep_resample.size());
        outfile.samples.front() = std::move(sweep_resample);
        outfile.save(qwqdsp_support::OutputFile("sweep_polyphase_resample_fir_down.wav"));
    }
    {
        qwqdsp_fx::PolyphaseUpsamplerFir resample;
        resample.Init(coeff, ratio);
        resample.Reset();
        auto sweep_resample = resample.Process(sweep);

        AudioFile<float> outfile;
        outfile.setNumChannels(1);
        outfile.setBitDepth(32);
        outfile.setSampleRate(infile.getSampleRate() * ratio);
        outfile.setNumSamplesPerChannel(sweep_resample.size());
        outfile.samples.front() = std::move(sweep_resample);
        outfile.save(qwqdsp_support::OutputFile("sweep_polyphase_resample_fir_up.wav"));
    }
}

// kPolyphaseCoeffs (N=512)
static constexpr std::array<float, 512> kPolyphaseCoeffs = {
    0.0000148924f,  -0.0000149728f, -0.0000122360f, -0.0000111857f, -0.0000103326f, -0.0000087969f, -0.0000061759f,
    -0.0000024478f, 0.0000020667f,  0.0000067966f,  0.0000110008f,  0.0000139155f,  0.0000148724f,  0.0000134568f,
    0.0000095923f,  0.0000036227f,  -0.0000037181f, -0.0000113583f, -0.0000180587f, -0.0000225736f, -0.0000238986f,
    -0.0000214247f, -0.0000151554f, -0.0000056478f, 0.0000057880f,  0.0000175926f,  0.0000277610f,  0.0000344565f,
    0.0000362615f,  0.0000323404f,  0.0000227281f,  0.0000084347f,  -0.0000086481f, -0.0000259991f, -0.0000408155f,
    -0.0000504395f, -0.0000528082f, -0.0000468521f, -0.0000327673f, -0.0000120976f, 0.0000123969f,  0.0000370850f,
    0.0000579833f,  0.0000713735f,  0.0000744409f,  0.0000657952f,  0.0000458501f,  0.0000168537f,  -0.0000172493f,
    -0.0000513990f, -0.0000801003f, -0.0000982749f, -0.0001021772f, -0.0000900384f, -0.0000625434f, -0.0000229079f,
    0.0000234360f,  0.0000695819f,  0.0001081195f,  0.0001323079f,  0.0001371954f,  0.0001205787f,  0.0000835360f,
    0.0000305100f,  -0.0000311887f, -0.0000923249f, -0.0001431112f, -0.0001747103f, -0.0001807559f, -0.0001585036f,
    -0.0001095683f, -0.0000399090f, 0.0000407806f,  0.0001204089f,  0.0001862469f,  0.0002269140f,  0.0002342971f,
    0.0002050612f,  0.0001414711f,  0.0000514161f,  -0.0000525046f, -0.0001546868f, -0.0002388288f, -0.0002904645f,
    -0.0002994003f, -0.0002615890f, -0.0001801620f, -0.0000653513f, 0.0000666925f,  0.0001961206f,  0.0003023086f,
    0.0003670962f,  0.0003778105f,  0.0003295999f,  0.0002266547f,  0.0000820720f,  -0.0000837226f, -0.0002457658f,
    -0.0003782916f, -0.0004587275f, -0.0004714784f, -0.0004107635f, -0.0002820895f, -0.0001019830f, 0.0001039924f,
    0.0003048127f,  0.0004685775f,  0.0005675199f,  0.0005825939f,  0.0005069712f,  0.0003477415f,  0.0001255497f,
    -0.0001279738f, -0.0003746002f, -0.0005752171f, -0.0006959213f, -0.0007136569f, -0.0006203740f, -0.0004250843f,
    -0.0001532897f, 0.0001561915f,  0.0004566739f,  0.0007005539f,  0.0008467648f,  0.0008675403f,  0.0007534582f,
    0.0005158015f,  0.0001858094f,  -0.0001892739f, -0.0005528352f, -0.0008473519f, -0.0010233644f, -0.0010476428f,
    -0.0009091642f, -0.0006219080f, -0.0002238289f, 0.0002279576f,  0.0006652466f,  0.0010189134f,  0.0012297223f,
    0.0012580589f,  0.0010910601f,  0.0007458448f,  0.0002682314f,  -0.0002731515f, -0.0007965515f, -0.0012193164f,
    -0.0014707810f, -0.0015038825f, -0.0013035816f, -0.0008906726f, -0.0003201208f, 0.0003259907f,  0.0009501005f,
    0.0014537263f,  0.0017528343f,  0.0017916051f,  0.0015524259f,  0.0010603193f,  0.0003809295f,  -0.0003879547f,
    -0.0011302464f, -0.0017288997f, -0.0020841375f, -0.0021297959f, -0.0018451229f, -0.0012600127f, -0.0004525580f,
    0.0004610233f,  0.0013428579f,  0.0020539672f,  0.0024759073f,  0.0025301307f,  0.0021919898f,  0.0014969376f,
    0.0005376425f,  -0.0005479495f, -0.0015961124f, -0.0024417235f, -0.0029439213f, -0.0030091219f, -0.0026076823f,
    -0.0017813548f, -0.0006399550f, 0.0006526954f,  0.0019018576f,  0.0029107883f,  0.0035112696f,  0.0035910557f,
    0.0031138749f,  0.0021285178f,  0.0007651446f,  -0.0007812195f, -0.0022780035f, -0.0034894887f, -0.0042132831f,
    -0.0043133261f, -0.0037441431f, -0.0025622135f, -0.0009220759f, 0.0009429299f,  0.0027530414f,  0.0042231938f,
    0.0051069686f,  0.0052367143f,  0.0045534898f,  0.0031217110f,  0.0011254980f,  -0.0011536293f, -0.0033752191f,
    -0.0051894317f, -0.0062906367f, -0.0064670521f, -0.0056386312f, -0.0038767767f, -0.0014018997f, 0.0014420021f,
    0.0042331459f,  0.0065322983f,  0.0079493108f,  0.0082061243f,  0.0071864982f,  0.0049641797f,  0.0018039838f,
    -0.0018660231f, -0.0055088421f, -0.0085528986f, -0.0104766381f, -0.0108914562f, -0.0096105935f, -0.0066929246f,
    -0.0024535194f, 0.0025628338f,  0.0076438999f,  0.0120016272f,  0.0148821729f,  0.0156800961f,  0.0140410616f,
    0.0099380493f,  0.0037088074f,  -0.0039531321f, -0.0120557537f, -0.0194096826f, -0.0247626184f, -0.0269516664f,
    -0.0250550801f, -0.0185241725f, -0.0072778772f, 0.0082535256f,  0.0271430901f,  0.0480408090f,  0.0692993630f,
    0.0891386083f,  0.1058303277f,  0.1178819344f,  0.1241975785f,  0.1241975785f,  0.1178819344f,  0.1058303277f,
    0.0891386083f,  0.0692993630f,  0.0480408090f,  0.0271430901f,  0.0082535256f,  -0.0072778772f, -0.0185241725f,
    -0.0250550801f, -0.0269516664f, -0.0247626184f, -0.0194096826f, -0.0120557537f, -0.0039531321f, 0.0037088074f,
    0.0099380493f,  0.0140410616f,  0.0156800961f,  0.0148821729f,  0.0120016272f,  0.0076438999f,  0.0025628338f,
    -0.0024535194f, -0.0066929246f, -0.0096105935f, -0.0108914562f, -0.0104766381f, -0.0085528986f, -0.0055088421f,
    -0.0018660231f, 0.0018039838f,  0.0049641797f,  0.0071864982f,  0.0082061243f,  0.0079493108f,  0.0065322983f,
    0.0042331459f,  0.0014420021f,  -0.0014018997f, -0.0038767767f, -0.0056386312f, -0.0064670521f, -0.0062906367f,
    -0.0051894317f, -0.0033752191f, -0.0011536293f, 0.0011254980f,  0.0031217110f,  0.0045534898f,  0.0052367143f,
    0.0051069686f,  0.0042231938f,  0.0027530414f,  0.0009429299f,  -0.0009220759f, -0.0025622135f, -0.0037441431f,
    -0.0043133261f, -0.0042132831f, -0.0034894887f, -0.0022780035f, -0.0007812195f, 0.0007651446f,  0.0021285178f,
    0.0031138749f,  0.0035910557f,  0.0035112696f,  0.0029107883f,  0.0019018576f,  0.0006526954f,  -0.0006399550f,
    -0.0017813548f, -0.0026076823f, -0.0030091219f, -0.0029439213f, -0.0024417235f, -0.0015961124f, -0.0005479495f,
    0.0005376425f,  0.0014969376f,  0.0021919898f,  0.0025301307f,  0.0024759073f,  0.0020539672f,  0.0013428579f,
    0.0004610233f,  -0.0004525580f, -0.0012600127f, -0.0018451229f, -0.0021297959f, -0.0020841375f, -0.0017288997f,
    -0.0011302464f, -0.0003879547f, 0.0003809295f,  0.0010603193f,  0.0015524259f,  0.0017916051f,  0.0017528343f,
    0.0014537263f,  0.0009501005f,  0.0003259907f,  -0.0003201208f, -0.0008906726f, -0.0013035816f, -0.0015038825f,
    -0.0014707810f, -0.0012193164f, -0.0007965515f, -0.0002731515f, 0.0002682314f,  0.0007458448f,  0.0010910601f,
    0.0012580589f,  0.0012297223f,  0.0010189134f,  0.0006652466f,  0.0002279576f,  -0.0002238289f, -0.0006219080f,
    -0.0009091642f, -0.0010476428f, -0.0010233644f, -0.0008473519f, -0.0005528352f, -0.0001892739f, 0.0001858094f,
    0.0005158015f,  0.0007534582f,  0.0008675403f,  0.0008467648f,  0.0007005539f,  0.0004566739f,  0.0001561915f,
    -0.0001532897f, -0.0004250843f, -0.0006203740f, -0.0007136569f, -0.0006959213f, -0.0005752171f, -0.0003746002f,
    -0.0001279738f, 0.0001255497f,  0.0003477415f,  0.0005069712f,  0.0005825939f,  0.0005675199f,  0.0004685775f,
    0.0003048127f,  0.0001039924f,  -0.0001019830f, -0.0002820895f, -0.0004107635f, -0.0004714784f, -0.0004587275f,
    -0.0003782916f, -0.0002457658f, -0.0000837226f, 0.0000820720f,  0.0002266547f,  0.0003295999f,  0.0003778105f,
    0.0003670962f,  0.0003023086f,  0.0001961206f,  0.0000666925f,  -0.0000653513f, -0.0001801620f, -0.0002615890f,
    -0.0002994003f, -0.0002904645f, -0.0002388288f, -0.0001546868f, -0.0000525046f, 0.0000514161f,  0.0001414711f,
    0.0002050612f,  0.0002342971f,  0.0002269140f,  0.0001862469f,  0.0001204089f,  0.0000407806f,  -0.0000399090f,
    -0.0001095683f, -0.0001585036f, -0.0001807559f, -0.0001747103f, -0.0001431112f, -0.0000923249f, -0.0000311887f,
    0.0000305100f,  0.0000835360f,  0.0001205787f,  0.0001371954f,  0.0001323079f,  0.0001081195f,  0.0000695819f,
    0.0000234360f,  -0.0000229079f, -0.0000625434f, -0.0000900384f, -0.0001021772f, -0.0000982749f, -0.0000801003f,
    -0.0000513990f, -0.0000172493f, 0.0000168537f,  0.0000458501f,  0.0000657952f,  0.0000744409f,  0.0000713735f,
    0.0000579833f,  0.0000370850f,  0.0000123969f,  -0.0000120976f, -0.0000327673f, -0.0000468521f, -0.0000528082f,
    -0.0000504395f, -0.0000408155f, -0.0000259991f, -0.0000086481f, 0.0000084347f,  0.0000227281f,  0.0000323404f,
    0.0000362615f,  0.0000344565f,  0.0000277610f,  0.0000175926f,  0.0000057880f,  -0.0000056478f, -0.0000151554f,
    -0.0000214247f, -0.0000238986f, -0.0000225736f, -0.0000180587f, -0.0000113583f, -0.0000037181f, 0.0000036227f,
    0.0000095923f,  0.0000134568f,  0.0000148724f,  0.0000139155f,  0.0000110008f,  0.0000067966f,  0.0000020667f,
    -0.0000024478f, -0.0000061759f, -0.0000087969f, -0.0000103326f, -0.0000111857f, -0.0000122360f, -0.0000149728f,
    0.0000148924f,
};

static void ResamplePolyphaseFir2() {
    constexpr int ratio = 8;

    AudioFile<float> infile;
    infile.load(qwqdsp_support::InputFile("sweep.wav"));
    auto& sweep = infile.samples.front();

    {
        qwqdsp_fx::PolyphaseDownsamplerFir resample;
        resample.Init(kPolyphaseCoeffs, ratio);
        resample.Reset();
        auto sweep_resample = resample.Process(sweep);

        AudioFile<float> outfile;
        outfile.setNumChannels(1);
        outfile.setBitDepth(32);
        outfile.setSampleRate(infile.getSampleRate() / ratio);
        outfile.setNumSamplesPerChannel(sweep_resample.size());
        outfile.samples.front() = std::move(sweep_resample);
        outfile.save(qwqdsp_support::OutputFile("sweep_polyphase_resample_fir_down2.wav"));
    }
    {
        qwqdsp_fx::PolyphaseUpsamplerFir resample;
        resample.Init(kPolyphaseCoeffs, ratio);
        resample.Reset();
        auto sweep_resample = resample.Process(sweep);

        AudioFile<float> outfile;
        outfile.setNumChannels(1);
        outfile.setBitDepth(32);
        outfile.setSampleRate(infile.getSampleRate() * ratio);
        outfile.setNumSamplesPerChannel(sweep_resample.size());
        outfile.samples.front() = std::move(sweep_resample);
        outfile.save(qwqdsp_support::OutputFile("sweep_polyphase_resample_fir_up2.wav"));
    }
}

int main() {
    ResampleIir();
    ResampleFir();
    ResamplePolyphaseFir();
    ResamplePolyphaseFir2();
}

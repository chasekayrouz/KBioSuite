# AlphaFold 3 in Kmatch

Kmatch predicts structures with Boltz-2. AlphaFold 3 (AF3) is not offered as a
predictor: it has not been tested in Kmatch, and you can't choose it in the app or
on the command line.

The **AlphaFold 3** card in **Dependencies & setup** is optional. Its **Install AF3
code** button (or `uv run kmatch setup af3-code` in a source build) installs AF3's
code only, for future testing: it downloads the pinned AF3 v3.0.4 source, creates its
Python environment and builds its chemical data. It downloads no model parameters and
no genetic databases. You don't need it to use Kmatch.

## Obtain parameters yourself

KBioSuite never includes or downloads AF3's model parameters (weights). Their terms
are separate from Kmatch and from AF3's source code.

The [official instructions](https://github.com/google-deepmind/alphafold3#obtaining-model-parameters)
provide a direct download from Google. Before downloading, read the
[parameter terms](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)
and the [prohibited-use policy](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_PROHIBITED_USE_POLICY.md)
and make sure you are authorized under them. The terms limit use to specified
non-commercial organizational activities and restrict redistribution. AF3's Apache
source license does not grant unrestricted use of the weights.

If you are authorized, run this in Terminal. It saves the parameters in `models/af3`
inside Kmatch's data folder (`~/Library/Application Support/KBioSuite/Kmatch/`). If
your Kmatch `config.json` sets `af3_models_override`, or you set `KMATCH_HOME`, use
that location instead.

```sh
mkdir -p "$HOME/Library/Application Support/KBioSuite/Kmatch/models/af3"
curl --fail --location --retry 5 --continue-at - \
  https://storage.googleapis.com/alphafold3/af3.bin.zst \
  --output "$HOME/Library/Application Support/KBioSuite/Kmatch/models/af3/af3.bin.zst.part"
mv "$HOME/Library/Application Support/KBioSuite/Kmatch/models/af3/af3.bin.zst.part" \
   "$HOME/Library/Application Support/KBioSuite/Kmatch/models/af3/af3.bin.zst"
```

Kmatch reads `af3.bin.zst` directly, so don't decompress it. If Google changes how
the parameters are obtained, follow the current official instructions; they list
`alphafold@google.com` for access questions. Don't use third-party copies of the
weights.

## Output terms

AF3 predictions are subject to the
[AF3 output terms](https://github.com/google-deepmind/alphafold3/blob/main/OUTPUT_TERMS_OF_USE.md).
Follow them when you export or share AF3 results.

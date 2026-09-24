# squonk2-train-test-split

![Data Manager Job: 2021.1](https://img.shields.io/badge/data%20manager%20job-2021.1-000000?labelColor=dc332e)

A Squonk2 Data Manager Job that splits a dataset into training, test and
validation subsets.

Splitting is either `random` or `scaffold`-based; scaffold splitting groups
molecules by their Murcko scaffold so that related structures do not straddle
the split, which gives a more honest estimate of how a model generalises.
Scaffolds can optionally be reduced by heavy-atom count (`hac`) or molecular
weight (`mw`).

The Job Definition lives in [`data-manager/jobs.yaml`](data-manager/jobs.yaml)
and belongs to the `im-virtual-screening` collection.

## The image

The image is built from [`Dockerfile`](Dockerfile) and published as
`informaticsmatters/train-test-split`. Dependencies are managed with Poetry;
`poetry.lock` is what the image installs, so a dependency change means
relocking with the Poetry version the `Dockerfile` pins.

Note that RDKit uses CalVer, so a caret constraint does not mean what it
usually does: `^2025.9.1` expands to `>=2025.9.1,<2026.0.0` and silently
excludes every 2026 release. Prefer a plain `>=` floor.

## Testing

Jobs are tested with [jote]:

```bash
poetry install --only dev
docker build -f Dockerfile . -t informaticsmatters/train-test-split:ci
jote --image-tag ci
```

That is exactly what CI does on every branch — it builds the image and runs the
real Job tests against it, rather than validating the Job Definition alone.
Nothing is pushed, so no registry credentials are needed.

## Releasing

1. Run the `publish-tag` workflow with the new tag.
2. Set the Job Definition `version` and `image.tag` to that same value.
3. Cut the matching Git tag.

Never reuse a container tag — the Data Manager caches any tag other than
`latest`/`stable` per Kubernetes node. See `docs/versioning.md` in the
[squonk2-jobs] umbrella repository.

## Architecture

Images are built for `linux/amd64`. The workflows this repository was forked
from declared `linux/amd64,linux/arm64`, but they never ran, and no arm64 image
has ever been published. Building RDKit, SciPy and scikit-learn for arm64 under
QEMU emulation is slow enough to want a dedicated runner rather than an
emulated one, so it is not enabled here.

[jote]: https://github.com/InformaticsMatters/squonk2-data-manager-job-tester
[squonk2-jobs]: https://github.com/InformaticsMatters/squonk2-jobs

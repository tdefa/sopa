from pathlib import Path


class WorkflowPaths:
    def __init__(self, sdata_path: str) -> None:
        self.sdata_path = Path(sdata_path)
        self.sdata_zgroup = self.sdata_path / ".zgroup"  # trick to fix snakemake ChildIOException
        self.raw = self.sdata_path.with_suffix(".qptiff")  # TODO: make it general

        self.polygons = self.sdata_path / "shapes" / "polygons"
        self.patches = self.sdata_path / "shapes" / "patches"
        self.table = self.sdata_path / "table"

        self.temp_dir = self.sdata_path.parent / f"{self.sdata_path.name}_temp"
        self.patches_dir = self.temp_dir / "patches"
        self.n_patches_path = self.sdata_path / ".n_patches"

        self.explorer_directory = self.sdata_path.with_suffix(".explorer")
        self.explorer_directory.mkdir(parents=True, exist_ok=True)

        self.explorer_experiment = self.explorer_directory / "experiment.xenium"

    def cells_paths(self, n: int):
        return [self.patches_dir / f"{i}.zarr.zip" for i in range(n)]


def _dump_arg(key: str, value):
    option = f"--{key.replace('_', '-')}"
    if isinstance(value, list):
        for v in value:
            yield from (option, str(v))
    elif value is True:
        yield option
    elif value is False:
        yield f"--no-{key.replace('_', '-')}"
    else:
        yield from (option, str(value))


def dump_args(args: dict) -> str:
    return " ".join((res for item in args.items() for res in _dump_arg(*item)))

from pathlib import Path
import subprocess
import shutil
class Wrapper:

    def __init__(self):
        self._path = "temporary_wrappers"

    def run(self):
        self.clone()
        self.copy_wrappers()
        self.teardown()

    def _get_path(self):
        return self._path
    repo_path = property(_get_path)

    def clone(self, tag=None):
        # clone the directory
        repo_path = Path(self.repo_path)

        if not repo_path.exists():
            print(f"Cloning sequana-wrappers into {repo_path}")
            repo_path.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run([
                "git", "clone", 
                "https://github.com/sequana/sequana-wrappers.git",
                "temporary_wrappers"
            ], check=True)
        if tag:
            print(f"Fetching tags and checking out tag: {tag}")
            subprocess.run(["git", "fetch", "--tags"], cwd=repo_path, check=True)
            subprocess.run(["git", "checkout", f"tags/{tag}", "-b", f"{tag}"], cwd=repo_path, check=True)
        else:
            print("Pulling latest changes from main") 
            subprocess.run(["git", "pull"], cwd=repo_path, check=True)

        return repo_path

    def copy_wrappers(self, destination="."):
        # copy the wrappers
        source_root = Path(self.repo_path)
        destination_root = Path(destination)
        print(f"Copying wrapper.py files to {destination_root}/wrappers")
        for wrapper_file in source_root.rglob("wrapper.py"):
            relative_path = wrapper_file.relative_to(source_root)
            dest_path = destination_root / relative_path
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(wrapper_file, dest_path)

    def teardown(self):
        # just suppress the clone
        shutil.rmtree(self.repo_path)


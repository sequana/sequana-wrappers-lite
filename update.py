from pathlib import Path
import subprocess
import shutil

import rich_click as click



class Wrapper:

    def __init__(self):
        self._path = "temporary_wrappers"

    def run(self, tag="main"):
        try:
            self.clone(tag=tag)
            self.copy_wrappers()
            self.teardown()
        except:
            self.teardown()

    def _get_path(self):
        return self._path
    repo_path = property(_get_path)

    def clone(self, tag="main"):
        # clone the directory
        repo_path = Path(self.repo_path)

        # main and latest version
        if not repo_path.exists():
            print(f"Cloning sequana-wrappers into {repo_path}")
            repo_path.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run([
                "git", "clone", 
                "https://github.com/sequana/sequana-wrappers.git",
                "temporary_wrappers"
            ], check=True)

        if tag and tag!= "main":
            print(f"Fetching tags and checking out tag: {tag}")
            subprocess.run(["git", "fetch", "--tags"], cwd=repo_path, check=True)
            subprocess.run(["git", "checkout", f"tags/{tag}", "-b", f"{tag}"], cwd=repo_path, check=True)
        else:
            print("Pulling latest changes from main")
            subprocess.run(["git", "pull"], cwd=repo_path, check=True)

        return repo_path

    def get_tags(self):
        result = subprocess.run(
            ["git", "tag"],
            cwd=".",
            capture_output=True,
            text=True,
            check=True
        )

        tags = result.stdout.strip().splitlines()
        return tags

    def get_origin_tags(self):
        # Get tags
        if not Path(self.repo_path).exists():
            self.clone()

        result = subprocess.run(
            ["git", "tag"],
            cwd=self.repo_path,
            capture_output=True,
            text=True,
            check=True
        )

        tags = result.stdout.strip().splitlines()
        # some tags are obviously wrong. we accept only those starting with vX.Y.Z
        tags = [x for x in tags if x.startswith('v')]
        return tags

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
        for wrapper_file in source_root.rglob("environment.yaml"):
            relative_path = wrapper_file.relative_to(source_root)
            dest_path = destination_root / relative_path
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(wrapper_file, dest_path)

    def teardown(self):
        # just suppress the clone
        shutil.rmtree(self.repo_path)




@click.command()
@click.option("--tag", type=click.STRING, required=True, help="must be set to 'main' or valid tag. type dummy one to get a valid list")
def main(**kwargs):

    w = Wrapper()
    w.clone()
    origin_tags = w.get_origin_tags()
    local_tags = w.get_tags()
    tag = kwargs['tag']

    if tag == "main":
        pass
    elif tag in local_tags:
        print(f"{tag} is already within the repository. Nothing to do.")    
        return
    elif tag not in origin_tags:
        print(f"must use a valid tag. choose one of {origin_tags}")
        return
    w.run(tag=tag)

    print("See README.md for the next steps")

    if tag == "main":
        help = """
git commit -m "update main branch"
git push origin main
"""
    else:
        help = f"""You should now do:

git commit -m "Add wrapper.py files for tag {tag} ."
git tag {tag}
git push origin main
git push origin {tag}
"""
    print(help)


if __name__ == '__main__':
    main()


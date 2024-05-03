# Silent Physics Animation Tasks for Learner Narration

This is the source code for silent physics animations for use in introductory physics classes.

> [!NOTE]
> - For easy viewing, the animations are available on [YouTube](https://www.youtube.com/channel/UC6jto5adixjNcKQ5PonID9w).
> - To download the animation video files, go to [Releases](https://github.com/SilentPhysicsAnimations/SilentPhysicsAnimations/releases) and click the desired zip-file.
> - To render the animations yourself, see [Rendering](#rendering).

The creation of these silent animation videos for physics education was
supported by the Icelandic Educational Materials Developement Fund ÞNS-2023-232672-1601.

## Customization

To customize any animation, one must edit its corresponding `.py` file and rerender (as specified in [Rendering](#rendering)).
Common attributes, such as colors and labels, are pulled to the top of the file where appropriate, for easy editing.
More advanced changes require editing the source code itself.

## Rendering

### Via Python

Before rendering the animations, ensure you have [Python](https://www.python.org/) installed (version 3.8 or later),
as well as [Manim Community](https://www.manim.community/) (v0.17.3). Then, clone this repository to your machine,
navigate to its root directory and run the following command:

```
python render.py
```

This will open a rendering wizard, where one can select the animations to render and the quality in which to render them.

The time to render will depend on hardware, but it is expected to take quite some time (upwards of 2 hours total if rendering in 4K).

### Via Docker

An alternative option to installing Python and Manim locally is using [Docker](https://www.docker.com/).
This will mitigate future version conflicts, and is the most reproducible option, although it might slow down rendering.
In this case the above command should instead be replaced with the following:

```
docker run --rm -it -v '.:/manim' manimcommunity/manim:v0.17.3 python render.py
```

This is otherwise identical to the alternative.

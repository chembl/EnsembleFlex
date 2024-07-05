# Original code comes from https://github.com/aidanjungo/StreamlitDirectoryPicker/blob/main/directorypicker.py
# Adapted with widget keys and a variation of session_state variable to allow for
# double usage (input and output) of the directory_picker in the same app.

import streamlit as st
from pathlib import Path


@st.experimental_fragment
def st_directory_picker_input(initial_path=Path(), key="key"):

    #st.markdown("#### Directory picker")

    if "inpath" not in st.session_state:
        st.session_state.inpath = initial_path.absolute()

    manual_input = st.text_input("Selected directory:", st.session_state.inpath, key=key)

    manual_input = Path(manual_input)
    if manual_input != st.session_state.inpath:
        st.session_state.inpath = manual_input
        st.rerun()

    _, col1, col2, col3, _ = st.columns([1, 1, 3, 1, 1])

    with col1:
        st.markdown("#")
        if st.button("⬅️", key=key+"left") and "inpath" in st.session_state:
            st.session_state.inpath = st.session_state.inpath.parent
            st.rerun()

    with col2:
        subdirectroies = [
            f.stem
            for f in st.session_state.inpath.iterdir()
            if f.is_dir()
            and (not f.stem.startswith(".") and not f.stem.startswith("__"))
        ]
        if subdirectroies:
            st.session_state.new_dir = st.selectbox(
                "Available subdirectories", sorted(subdirectroies), key=key+"subdirectory"
            )
        else:
            st.markdown("#")
            st.markdown(
                "<font color='#FF0000'>No subdir</font>", unsafe_allow_html=True
            )

    with col3:
        if subdirectroies:
            st.markdown("#")
            if st.button("➡️", key=key+"right") and "inpath" in st.session_state:
                st.session_state.inpath = Path(
                    st.session_state.inpath, st.session_state.new_dir
                )
                st.rerun()

    return st.session_state.inpath


@st.experimental_fragment
def st_directory_picker_output(initial_path=Path(), key="key"):

    #st.markdown("#### Directory picker")

    if "outpath" not in st.session_state:
        st.session_state.outpath = initial_path.absolute()

    manual_input = st.text_input("Selected directory:", st.session_state.outpath, key=key)

    manual_input = Path(manual_input)
    if manual_input != st.session_state.outpath:
        st.session_state.outpath = manual_input
        st.rerun()

    _, col1, col2, col3, _ = st.columns([1, 1, 3, 1, 1])

    with col1:
        st.markdown("#")
        if st.button("⬅️", key=key+"left") and "outpath" in st.session_state:
            st.session_state.outpath = st.session_state.outpath.parent
            st.rerun()

    with col2:
        subdirectroies = [
            f.stem
            for f in st.session_state.outpath.iterdir()
            if f.is_dir()
            and (not f.stem.startswith(".") and not f.stem.startswith("__"))
        ]
        if subdirectroies:
            st.session_state.new_dir = st.selectbox(
                "Available subdirectories", sorted(subdirectroies), key=key+"subdirectory"
            )
        else:
            st.markdown("#")
            st.markdown(
                "<font color='#FF0000'>No subdir</font>", unsafe_allow_html=True
            )

    with col3:
        if subdirectroies:
            st.markdown("#")
            if st.button("➡️", key=key+"right") and "outpath" in st.session_state:
                st.session_state.outpath = Path(
                    st.session_state.outpath, st.session_state.new_dir
                )
                st.rerun()

    return st.session_state.outpath


@st.experimental_fragment
def st_directory_picker_input_liganded(initial_path=Path(), key="key"):

    #st.markdown("#### Directory picker")

    if "inpath_liganded" not in st.session_state:
        st.session_state.inpath_liganded = initial_path.absolute()

    manual_input = st.text_input("Selected directory:", st.session_state.inpath_liganded, key=key)

    manual_input = Path(manual_input)
    if manual_input != st.session_state.inpath_liganded:
        st.session_state.inpath_liganded = manual_input
        st.rerun()

    _, col1, col2, col3, _ = st.columns([1, 1, 3, 1, 1])

    with col1:
        st.markdown("#")
        if st.button("⬅️", key=key+"left") and "inpath" in st.session_state:
            st.session_state.inpath_liganded = st.session_state.inpath_liganded.parent
            st.rerun()

    with col2:
        subdirectroies = [
            f.stem
            for f in st.session_state.inpath_liganded.iterdir()
            if f.is_dir()
            and (not f.stem.startswith(".") and not f.stem.startswith("__"))
        ]
        if subdirectroies:
            st.session_state.new_dir = st.selectbox(
                "Available subdirectories", sorted(subdirectroies), key=key+"subdirectory"
            )
        else:
            st.markdown("#")
            st.markdown(
                "<font color='#FF0000'>No subdir</font>", unsafe_allow_html=True
            )

    with col3:
        if subdirectroies:
            st.markdown("#")
            if st.button("➡️", key=key+"right") and "inpath" in st.session_state:
                st.session_state.inpath_liganded = Path(
                    st.session_state.inpath_liganded, st.session_state.new_dir
                )
                st.rerun()

    return st.session_state.inpath_liganded

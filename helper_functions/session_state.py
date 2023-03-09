import streamlit as st


class States():
    def initialise_state(self, state_dict):
        for k,v in state_dict.items():
            if k not in st.session_state:
                st.session_state[k] = v

    def save_state(self, state_dict):
        for k,v in state_dict.items():
            st.session_state[k] = v
        return

ss = States()
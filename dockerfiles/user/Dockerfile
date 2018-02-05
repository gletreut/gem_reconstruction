FROM root/gem:16.04

# add new sudo user
USER 0
ENV USERNAME yourusername
ENV UID youruid
ENV GROUPNAME yourgroupname
ENV GID yourgid
RUN groupadd --gid $GID $GROUPNAME && \
    useradd -g $GROUPNAME -u $UID -m $USERNAME && \
    echo "$USERNAME:$USERNAME" | chpasswd && \
    usermod --shell /bin/bash $USERNAME && \
    usermod -aG sudo $USERNAME && \
    echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/$USERNAME && \
    chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME
WORKDIR /home/$USERNAME

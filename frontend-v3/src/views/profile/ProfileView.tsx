import { Button, Card, Elevation } from "@blueprintjs/core";
import { useProfileStore } from "../../modules/profile";
import styles from "./Profile.module.scss";
import { EditProfileDialog } from "./EditProfileDialog";
import { useState } from "react";
import { EditPasswordDialog } from "./EditPasswordDialog";

export default function ProfileView() {
  const userProfile = useProfileStore((state) => state.userProfile);
  const [editProfileDialogOpen, setEditProfileDialogOpen] = useState<boolean>(false);
  const [editPasswordDialogOpen, setEditPasswordDialogOpen] = useState<boolean>(false);

  return (
    <div className={styles.container}>
      <Card interactive={false} elevation={Elevation.TWO}>
        <h2>User Profile</h2>
        <h3>Name:</h3>
        <p>{userProfile?.name}</p>
        <h3>Email:</h3>
        <p>{userProfile?.email}</p>
        <div className={styles.actions}>
          <Button text="Edit profile" onClick={() => setEditProfileDialogOpen(true)} />
          <Button text="Change password" onClick={() => setEditPasswordDialogOpen(true)} />
        </div>
      </Card>
      <EditProfileDialog isOpen={editProfileDialogOpen} handleClose={() => setEditProfileDialogOpen(false)} />
      <EditPasswordDialog isOpen={editPasswordDialogOpen} handleClose={() => setEditPasswordDialogOpen(false)} />
    </div>
  );
}

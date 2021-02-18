import styles from "./GroupsListView.module.scss";
import shallow from "zustand/shallow";
import { useEffect, useState } from "react";
import { Alert, Button, Card, Classes, Elevation, Icon, InputGroup, Intent, Tag } from "@blueprintjs/core";
import { throttle } from "lodash-es";
import { AddGroupDialog } from "./AddGroupDialog";
import { EditGroupDialog } from "./EditGroupDialog";
import { useGroupsStore } from "modules/groups";
import { IGroup } from "modules/groups/models";
import Masonry from "react-masonry-css";
import "./masonry.scss";
import { useProfileStore } from "../../modules/profile";
import { useHistory } from "react-router-dom";

export function GroupsListView() {
  const history = useHistory();
  const userProfile = useProfileStore((state) => state.userProfile);
  const { ids, entities, deleteGroup, joinGroup, groupsTags } = useGroupsStore(
    (state) => ({
      ids: state.ids,
      entities: state.entities,
      deleteGroup: state.deleteGroup,
      joinGroup: state.joinGroup,
      groupsTags: state.groupsTags,
    }),
    shallow
  );
  const [filterValue, setFilterValue] = useState<string>("");
  const [editDialogOpen, setEditDialogOpen] = useState<boolean>(false);
  const [addDialogOpen, setAddDialogOpen] = useState<boolean>(false);
  const [alertOpen, setAlertOpen] = useState<boolean>(false);
  const [activeItem, setActiveItem] = useState<IGroup | null>(null);

  const groups = ids.map((id) => entities[id]);
  const data = groups.filter((item) => {
    const filter = filterValue.toLowerCase();
    return (
      item.name.toLowerCase().includes(filter) || (item.description && item.description.toLowerCase().includes(filter))
    );
  });

  const handleFilterChange = (event: React.FormEvent<HTMLElement>) =>
    setFilterValue((event.target as HTMLInputElement).value);

  const editAction = (item: IGroup) => {
    setActiveItem(item);
    setEditDialogOpen(true);
  };

  const deleteAction = async () => {
    if (activeItem) {
      await deleteGroup(activeItem.id);
    }
    setAlertOpen(false);
  };

  return (
    <div className={styles.container}>
      <span className={styles.toolbar}>
        <h2>Groups</h2>
        <Button text="Add Group" icon="add" intent={Intent.PRIMARY} onClick={() => setAddDialogOpen(true)} />
      </span>
      <InputGroup
        asyncControl={true}
        leftIcon="filter"
        rightElement={
          filterValue ? <Button icon="cross" minimal={true} onClick={() => setFilterValue("")} /> : undefined
        }
        onChange={throttle(handleFilterChange, 200, { leading: false })}
        placeholder="Filter groups..."
        value={filterValue}
        className={styles.filter}
      />
      <Masonry
        breakpointCols={{ default: 4, 1000: 3, 700: 2, 400: 1 }}
        className="masonry-grid"
        columnClassName="masonry-grid_column"
      >
        {data.map((group) => {
          const userIds = group.members.filter((item) => item.is_active).map((member) => member.user_id);
          const isMember = userIds.includes(userProfile!.id);
          return (
            <Card key={group.id} interactive={false} elevation={Elevation.TWO}>
              <h3>{group.name}</h3>
              {group.url && (
                <span className="bp3-text-small">
                  <Icon icon="globe-network" />{" "}
                  <a href={group.url} target="_blank" rel="noreferrer">
                    {group.url}
                  </a>
                </span>
              )}
              {group.description && <h5>{group.description}</h5>}
              {group.tags && (
                <p>
                  {group.tags.map((tag) => {
                    return (
                      <Tag key={tag} minimal={true} className={styles.tag}>
                        {tag}
                      </Tag>
                    );
                  })}
                </p>
              )}
              <div className={Classes.DIALOG_FOOTER_ACTIONS}>
                {(isMember || userProfile!.is_admin) && (
                  <Button
                    text="Open"
                    intent={Intent.PRIMARY}
                    onClick={() => history.push(`/${group.id}`)}
                  />
                )}
                {!isMember && group.is_open && <Button text="Join" onClick={() => joinGroup(group.id)} />}
                {isMember && userProfile!.is_admin && <Button text="Edit" onClick={() => editAction(group)} />}
                {isMember && userProfile!.is_admin && (
                  <Button
                    text="Delete"
                    intent={Intent.DANGER}
                    onClick={() => {
                      setActiveItem(group);
                      setAlertOpen(true);
                    }}
                  />
                )}
              </div>
            </Card>
          );
        })}
      </Masonry>
      {activeItem && (
        <EditGroupDialog
          group={activeItem}
          groupsTags={groupsTags!}
          isOpen={editDialogOpen}
          handleClose={() => {
            setEditDialogOpen(false);
            setActiveItem(null);
          }}
        />
      )}
      <AddGroupDialog groupsTags={groupsTags!} isOpen={addDialogOpen} handleClose={() => setAddDialogOpen(false)} />
      <Alert
        className={Classes.DARK}
        cancelButtonText="Cancel"
        canEscapeKeyCancel={true}
        confirmButtonText="Delete"
        icon="trash"
        intent={Intent.DANGER}
        isOpen={alertOpen}
        onCancel={() => setAlertOpen(false)}
        onConfirm={deleteAction}
      >
        <p>Are you sure you want to delete the group?</p>
      </Alert>
    </div>
  );
}

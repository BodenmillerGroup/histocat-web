import { Button, Checkbox, Classes, Dialog, FormGroup, InputGroup, Intent, MenuItem } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { useGroupsStore } from "modules/groups";
import { IGroupCreate } from "modules/groups/models";
import { ItemPredicate, ItemRenderer, MultiSelect } from "@blueprintjs/select";
import shallow from "zustand/shallow";
import { useState } from "react";

type AddGroupDialogProps = {
  isOpen: boolean;
  handleClose(): void;
  groupsTags: string[];
};

export function AddGroupDialog(props: AddGroupDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const { createGroup } = useGroupsStore(
    (state) => ({
      createGroup: state.createGroup,
    }),
    shallow
  );
  const [items, setItems] = useState<string[]>(props.groupsTags);
  const [tags, setTags] = useState<string[]>([]);
  const [createdItems, setCreatedItems] = useState<string[]>([]);

  const onSubmit = async (values: any) => {
    const params: IGroupCreate = {
      name: values.name,
      description: values.description,
      url: values.url,
      is_open: values.isPublic,
      tags: tags,
    };
    await createGroup(params);
    setTags([]);
    props.handleClose();
  };

  const renderCreateTagOption = (query: string, active: boolean, handleClick: React.MouseEventHandler<HTMLElement>) => (
    <MenuItem
      icon="add"
      text={`Create "${query}"`}
      active={active}
      onClick={handleClick}
      shouldDismissPopover={false}
    />
  );

  function createTag(title: string): string {
    return title;
  }

  function areTagsEqual(a: string, b: string) {
    // Compare only the titles (ignoring case) just for simplicity.
    return a.toLowerCase() === b.toLowerCase();
  }

  function arrayContainsTag(tags: string[], tagToFind: string): boolean {
    return tags.some((tag: string) => tag === tagToFind);
  }

  function addTagToArray(tags: string[], tagToAdd: string) {
    return [...tags, tagToAdd];
  }

  function deleteTagFromArray(tags: string[], tagToDelete: string) {
    return tags.filter((tag) => tag !== tagToDelete);
  }

  function maybeAddCreatedTagToArrays(
    items: string[],
    createdItems: string[],
    tag: string
  ): { createdItems: string[]; items: string[] } {
    const isNewlyCreatedItem = !arrayContainsTag(items, tag);
    return {
      createdItems: isNewlyCreatedItem ? addTagToArray(createdItems, tag) : createdItems,
      // Add a created film to `items` so that the film can be deselected.
      items: isNewlyCreatedItem ? addTagToArray(items, tag) : items,
    };
  }

  function maybeDeleteCreatedTagFromArrays(
    items: string[],
    createdItems: string[],
    tag: string
  ): { createdItems: string[]; items: string[] } {
    const wasItemCreatedByUser = arrayContainsTag(createdItems, tag);

    // Delete the item if the user manually created it.
    return {
      createdItems: wasItemCreatedByUser ? deleteTagFromArray(createdItems, tag) : createdItems,
      items: wasItemCreatedByUser ? deleteTagFromArray(items, tag) : items,
    };
  }

  const renderTag = (tag: string) => tag;

  // NOTE: not using Films.itemRenderer here so we can set icons.
  const renderFilm: ItemRenderer<string> = (film, { modifiers, handleClick }) => {
    if (!modifiers.matchesPredicate) {
      return null;
    }
    return (
      <MenuItem
        active={modifiers.active}
        icon={isTagSelected(film) ? "tick" : "blank"}
        key={film}
        onClick={handleClick}
        text={film}
        shouldDismissPopover={false}
      />
    );
  };

  const handleTagRemove = (_tag: React.ReactNode, index: number) => {
    deselectTag(index);
  };

  const getSelectedTagIndex = (tag: string) => {
    return tags.indexOf(tag);
  };

  const isTagSelected = (tag: string) => {
    return getSelectedTagIndex(tag) !== -1;
  };

  const selectTag = (tag: string) => {
    selectTags([tag]);
  };

  const selectTags = (tagsToSelect: string[]) => {
    let nextCreatedItems = createdItems.slice();
    let nextTags = tags.slice();
    let nextItems = items.slice();

    tagsToSelect.forEach((tag) => {
      const results = maybeAddCreatedTagToArrays(nextItems, nextCreatedItems, tag);
      nextItems = results.items;
      nextCreatedItems = results.createdItems;
      // Avoid re-creating an item that is already selected (the "Create
      // Item" option will be shown even if it matches an already selected
      // item).
      nextTags = !arrayContainsTag(nextTags, tag) ? [...nextTags, tag] : nextTags;
    });

    setCreatedItems(nextCreatedItems);
    setTags(nextTags);
    setItems(nextItems);
  };

  const deselectTag = (index: number) => {
    const tag = tags[index];
    const { createdItems: nextCreatedItems, items: nextItems } = maybeDeleteCreatedTagFromArrays(
      items,
      createdItems,
      tag
    );

    // Delete the item if the user manually created it.
    setCreatedItems(nextCreatedItems);
    setTags(tags.filter((_tag, i) => i !== index));
    setItems(nextItems);
  };

  const handleTagSelect = (tag: string) => {
    if (!isTagSelected(tag)) {
      selectTag(tag);
    } else {
      deselectTag(getSelectedTagIndex(tag));
    }
  };

  const handleTagsPaste = (tags: string[]) => {
    // On paste, don't bother with deselecting already selected values, just
    // add the new ones.
    selectTags(tags);
  };

  const filterTag: ItemPredicate<string> = (query, tag, _index, exactMatch) => {
    const normalizedTitle = tag.toLowerCase();
    const normalizedQuery = query.toLowerCase();

    if (exactMatch) {
      return normalizedTitle === normalizedQuery;
    } else {
      return tag.indexOf(normalizedQuery) >= 0;
    }
  };

  const handleClear = () => setTags([]);
  const clearButton = tags.length > 0 ? <Button icon="cross" minimal={true} onClick={handleClear} /> : undefined;

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Add Group"
      usePortal={true}
      isOpen={props.isOpen}
      className="bp3-dark"
      canOutsideClickClose={false}
    >
      <form onSubmit={handleSubmit(onSubmit)}>
        <div className={Classes.DIALOG_BODY}>
          <FormGroup
            label="Name"
            labelFor="name-input"
            labelInfo="(required)"
            intent="danger"
            helperText={errors?.name?.message}
          >
            <InputGroup
              id="name-input"
              name="name"
              placeholder="Enter group name"
              inputRef={register({
                required: "Name is required",
              })}
            />
          </FormGroup>

          <FormGroup label="Description" labelFor="description-input">
            <InputGroup
              id="description-input"
              name="description"
              placeholder="Enter group description"
              inputRef={register({})}
            />
          </FormGroup>

          <FormGroup label="URL" labelFor="url-input">
            <InputGroup id="url-input" name="url" placeholder="Enter group URL" inputRef={register({})} />
          </FormGroup>

          <FormGroup label="Tags">
            <MultiSelect
              createNewItemFromQuery={createTag}
              createNewItemRenderer={renderCreateTagOption}
              itemRenderer={renderFilm}
              itemPredicate={filterTag}
              itemsEqual={areTagsEqual}
              // we may customize the default filmSelectProps.items by
              // adding newly created items to the list, so pass our own
              items={items}
              noResults={<MenuItem disabled={true} text="No results." />}
              onItemSelect={handleTagSelect}
              onItemsPaste={handleTagsPaste}
              popoverProps={{ minimal: true }}
              tagRenderer={renderTag}
              tagInputProps={{
                onRemove: handleTagRemove,
                rightElement: clearButton,
                tagProps: { minimal: true },
              }}
              selectedItems={tags}
              resetOnSelect={true}
              fill={true}
            />
          </FormGroup>

          <Checkbox name="isPublic" label="Public group" inline={true} defaultChecked={false} inputRef={register({})} />
        </div>
        <div className={Classes.DIALOG_FOOTER}>
          <div className={Classes.DIALOG_FOOTER_ACTIONS}>
            <Button onClick={props.handleClose} text="Cancel" />
            <Button type="reset" text="Reset" onClick={() => setTags([])} />
            <Button type="submit" intent={Intent.PRIMARY} text="Save" />
          </div>
        </div>
      </form>
    </Dialog>
  );
}

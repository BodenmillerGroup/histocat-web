import { Button, Checkbox, Classes, Dialog, FormGroup, InputGroup, Intent, MenuItem } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { useGroupsStore } from "modules/groups";
import { IGroup, IGroupUpdate } from "modules/groups/models";
import { useState } from "react";
import { ItemPredicate, ItemRenderer, MultiSelect } from "@blueprintjs/select";

type EditModelDialogProps = {
  isOpen: boolean;
  handleClose(): void;
  group: IGroup;
  groupsTags: string[];
};

export function EditGroupDialog(props: EditModelDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const updateGroup = useGroupsStore((state) => state.updateGroup);
  const [items, setItems] = useState<string[]>(props.groupsTags);
  const [tags, setTags] = useState<string[]>(props.group.tags);
  const [createdItems, setCreatedItems] = useState<string[]>([]);

  const onSubmit = async (values: any) => {
    // form is valid
    const params: IGroupUpdate = {
      name: values.name,
      description: values.description,
      url: values.url,
      is_open: values.isOpen,
      tags: tags,
    };
    await updateGroup(props.group.id, params);
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
      // Add a created tag to `items` so that the tag can be deselected.
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

  const renderItem: ItemRenderer<string> = (tag, { modifiers, handleClick }) => {
    if (!modifiers.matchesPredicate) {
      return null;
    }
    return (
      <MenuItem
        active={modifiers.active}
        icon={isTagSelected(tag) ? "tick" : "blank"}
        key={tag}
        onClick={handleClick}
        text={tag}
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
      title="Edit Group"
      usePortal={true}
      isOpen={props.isOpen}
      className={Classes.DARK}
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
              defaultValue={props.group.name}
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
              defaultValue={props.group.description}
              inputRef={register({})}
            />
          </FormGroup>

          <FormGroup label="URL" labelFor="url-input">
            <InputGroup
              id="url-input"
              name="url"
              placeholder="Enter group URL"
              defaultValue={props.group.url}
              inputRef={register({})}
            />
          </FormGroup>

          <FormGroup label="Tags">
            <MultiSelect
              createNewItemFromQuery={createTag}
              createNewItemRenderer={renderCreateTagOption}
              itemRenderer={renderItem}
              itemPredicate={filterTag}
              itemsEqual={areTagsEqual}
              // we may customize the default tagSelectProps.items by
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

          <Checkbox
            name="isOpen"
            label="Public group"
            inline={true}
            defaultChecked={props.group.is_open}
            inputRef={register({})}
          />
        </div>
        <div className={Classes.DIALOG_FOOTER}>
          <div className={Classes.DIALOG_FOOTER_ACTIONS}>
            <Button onClick={props.handleClose} text="Cancel" />
            <Button type="reset" text="Reset" onClick={() => setTags(props.group.tags)} />
            <Button type="submit" intent={Intent.PRIMARY} text="Save" />
          </div>
        </div>
      </form>
    </Dialog>
  );
}

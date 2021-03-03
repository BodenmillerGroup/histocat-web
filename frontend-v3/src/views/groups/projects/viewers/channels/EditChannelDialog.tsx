import { Button, Classes, Dialog, FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { useForm } from "react-hook-form";
import { IChannel, IChannelUpdate } from "modules/projects/models";
import { useProjectsStore } from "modules/projects";
import shallow from "zustand/shallow";

type EditChannelDialogProps = {
  isOpen: boolean;
  handleClose(): void;
  channel: IChannel;
};

export function EditChannelDialog(props: EditChannelDialogProps) {
  const { register, errors, handleSubmit } = useForm();
  const { activeAcquisitionId, updateChannel } = useProjectsStore(
    (state) => ({
      activeAcquisitionId: state.activeAcquisitionId,
      updateChannel: state.updateChannel,
    }),
    shallow
  );

  const onSubmit = async (values: any) => {
    // form is valid
    const params: IChannelUpdate = {
      name: props.channel.name,
      customLabel: values.label,
    };
    await updateChannel(activeAcquisitionId!, params);
    props.handleClose();
  };

  return (
    <Dialog
      icon="edit"
      onClose={props.handleClose}
      title="Edit Channel"
      usePortal={true}
      isOpen={props.isOpen}
      className={Classes.DARK}
      canOutsideClickClose={false}
    >
      <form onSubmit={handleSubmit(onSubmit)}>
        <div className={Classes.DIALOG_BODY}>
          <FormGroup label="Label" labelFor="label-input" labelInfo="(required)">
            <InputGroup
              id="label-input"
              name="label"
              placeholder="Enter channel label"
              defaultValue={props.channel.customLabel}
              inputRef={register({ required: "label is required" })}
            />
          </FormGroup>
        </div>
        <div className={Classes.DIALOG_FOOTER}>
          <div className={Classes.DIALOG_FOOTER_ACTIONS}>
            <Button onClick={props.handleClose} text="Cancel" />
            <Button type="reset" text="Reset" />
            <Button type="submit" intent={Intent.PRIMARY} text="Save" />
          </div>
        </div>
      </form>
    </Dialog>
  );
}
